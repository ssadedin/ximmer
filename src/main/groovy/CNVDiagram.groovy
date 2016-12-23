import java.awt.BasicStroke;
import java.text.NumberFormat

import jsr166y.ForkJoinPool;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction

import graxxia.Matrix
import graxxia.Stats
import groovy.transform.CompileStatic;
import groovy.util.logging.Log;
import groovy.xml.Namespace;
import groovyx.gpars.GParsPool;
import htsjdk.samtools.SAMFileReader

/**
 * Draws a diagram of the region of each CNV in a provided CNV summary report
 * 
 * @author simon
 *
 */
@Log
class CNVDiagram {
    
    def err = System.err
    
    final static NumberFormat numberFormat = NumberFormat.numberInstance
    
    static {
        numberFormat.maximumFractionDigits=1
        numberFormat.minimumFractionDigits=0
    }
    
    Map<String,SampleInfo> sampleInfo
    
    Regions cnvs
    
    /**
     * The samples for which CNVs will be plotted
     */
    List<String> samples 
    
    /**
     * All the samples that we have available
     */
    List<String> allSamples 
    
    /**
     * Coverage information from BEDTools
     */
    Map<String,SAM> bams
    
    /**
     * Calls from each invididual CNV caller
     */
    Map<String, RangedData> cnvCalls
    
    Matrix readCounts = null
    
    Matrix normCounts = null
    
    Map<String,Double> means = null
    
    RefGenes geneDb = null
    
    Regions targetRegions = null
    
    Map<String,Regions> vcfs = Collections.synchronizedMap([:])
    
    /**
     * Amplicons covering region (for amplicon based sequencing)
     */
    Regions amplicons = null
    
    /**
     * Amplicon read counts (for amplicon based sequencing)
     */
    Matrix ampliconCounts = null
    
    int concurrency = 1
    
    boolean extendedMeanEstimator = false
    
    boolean liteMeanEstimator = false
    
    boolean noAmpliconLabels = false
    
    String gatkMeanEstimator = null
    
    /**
     * Trying to display too many target regions in plots makes them unreadable.
     * This count limits the number of flanking regions that will be drawn next
     * to each deletion in cases where the affected gene is very large or spans a
     * huge number of target regions
     */
    int maxTargetCount = 20
    
    CNVDiagram(Regions cnvs, Map<String,SampleInfo> sampleInfo, Map<String, RangedData> cnvCalls, Regions targetRegions, List<String> samples=null) {
        
        this.cnvs = cnvs;
        
        this.cnvs.each { cnv ->
            if(cnv.sample == null)
                cnv.sample = cnv.extra
        }
        
        this.sampleInfo = sampleInfo
        this.samples = samples  == null ? this.cnvs*.sample : samples
        
        allSamples = (sampleInfo.keySet() as List).grep { 
            sampleInfo[it].files.bam
        }
        this.cnvCalls = cnvCalls
        this.targetRegions = targetRegions 
    }
    
    void loadBAMs() {
        List<SAM> sampleBams = allSamples.collect { String sample ->
            
            if(!sampleInfo[sample].files.bam)
                throw new IllegalArgumentException("Sample $sample does not have a BAM file provided")
                
            new SAM(sampleInfo[sample].files.bam[0]) 
        }
        
        bams = Collections.synchronizedMap([allSamples, sampleBams].transpose().collectEntries())
    }
    
    void loadMeans() {
        
        if(!bams)
            loadBAMs()
        
        GParsPool.withPool(concurrency) {
            log.info "Calculating coverage means for: $allSamples ..."
            
            Regions meanRegions = targetRegions
            
            if(liteMeanEstimator) {
                // Search for the smallest chromosome that we can that is in the target regions
                meanRegions = findLiteMeanRegions()
            }
            
            if(gatkMeanEstimator) {
                GATKMeanEstimator estimator = new GATKMeanEstimator(this.gatkMeanEstimator)
                log.info "Discovered ${estimator.intervalFiles.size()} samples with GATK coverage information"
                means = estimator.calculateMeans(allSamples)
            }
            else
            if(extendedMeanEstimator || liteMeanEstimator) {
                if(extendedMeanEstimator)
                    log.info "Using extended mean coverage estimation ..."
                else
                    log.info "Using lite mean coverage estimation ..."
                means = Collections.synchronizedMap([allSamples, 
                                                     allSamples.collectParallel { bams[it].coverageStatistics(meanRegions).mean }]
                                                               .transpose().collectEntries())
            }
            else {
                int targetSize = targetRegions.size()
                means = Collections.synchronizedMap([allSamples, allSamples.collectParallel { 
                        log.info "Counting reads for $it"; 
                        if(!bams[it])
                            System.err.println "No BAM file provided for $it"
                            
                        long count=0; bams[it].eachRecord { count += it.readLength }; 
                        count/(double)targetSize 
                    }
                ].transpose().collectEntries())
            }
        }
        log.info "Sample means are $means"
    }
    
    void draw(String outputFileBase, int width, int height) {
        
        if(!cnvs.count { 1 }) // if no CNV calls, short circuit loading all the coverage info
            return
        
        if(!bams)
            loadBAMs()
        
        if(!means) 
            loadMeans()
        
        def callers = cnvCalls.keySet() as List
        List<String> palette = ["red","green","orange","blue","gray","magenta","yellow","cyan","black"]
        def colors = [ callers, palette[0..<callers.size()] ].transpose().collectEntries()
        
        Closure processCnv = { cnv ->
            Utils.time("Draw CNV $cnv.chr:$cnv.from-$cnv.to") {
                drawCNV(cnv, outputFileBase, width, height, callers, colors)
            }
        }

        if(concurrency > 1) {
            GParsPool.withPool(concurrency) {
                cnvs.eachParallel(processCnv)
            }
        }
        else {
            cnvs.each(processCnv)
        }
    }

    void drawCNV(Region cnv, String outputFileBase, int width, int height, List callers, Map colors) {
        
        if(this.samples && !(cnv.sample in this.samples)) {
            println "Skip CNV $cnv.chr:$cnv.from-$cnv.to for sample $cnv.sample"
            return
        }
        
        println "Call $cnv.chr:$cnv.from-$cnv.to for sample $cnv.sample"

        String imageFileName = outputFileBase.replaceAll('.png$','') + "_${cnv.chr}_${cnv.from}_${cnv.to}_${cnv.sample}.png"
        File imageFile = new File(imageFileName)
        if(imageFile.exists() && imageFile.length() > 0) {
            println "Skip $imageFileName because file already exists"
            return
        }
        
        SAM bam = bams[cnv.sample]

        List<String> genes = geneDb.getGenes(cnv).grep { gene -> geneDb.getExons(gene).overlaps(this.targetRegions) }

        println "Genes are $genes"
        List targets = null
        
        // If one or more genes is overlapped, display the whole gene(s)
        if(!genes.empty) {
    
            Regions exons = genes.sum { geneDb.getExons(it) }.reduce().grep { it.chr == cnv.chr } as Regions
    
            println "Overlapping exons are " + exons + " from " + exons[0].from + " to " + exons[-1].to
    
            println "Exons are " + exons.collect { it.from + "-" + it.to }.join(",")
            
            // Find the overlapping target regions from the coverage file
            targets = targetRegions.getOverlaps(cnv.chr, exons[0].from, exons[-1].to)
            println "Overlapping targets are " + targets.collect { it.from+"-"+it.to }.join(",")
            
            // If there are too many targets, clip at a fixed number upstream and downstream
            // from the actual CNV
            if(targets.size()>maxTargetCount) {
                int cnvIndex = targets.findIndexOf {GRange.overlaps(it, cnv.range)}
                int lastCnvIndex = targets.findLastIndexOf {GRange.overlaps(it, cnv.range)}
                int flankingTargets = (int)Math.max(2,maxTargetCount - (lastCnvIndex - cnvIndex))/2
                targets = targets[(Math.max(0,cnvIndex-flankingTargets)..Math.min(lastCnvIndex+flankingTargets, targets.size()-1))]
            }
        }
        else { // no gene overlapped, just show the region of the CNV itself
            targets = targetRegions.getOverlaps(cnv.chr, cnv.from, cnv.to)
        }
        
        def froms = targets*.from
        def tos = targets*.to
        
        if(froms.size()==0 || tos.size()==0) {
            println  "ERROR: no targets overlapped by $cnv.chr:$cnv.from-$cnv.to"
            return
        }
        
        def displayRegion = new Region(cnv.chr, froms[0]..tos[-1])
        
        double yMax = 2.0d

        graxxia.Drawing d = new graxxia.Drawing(
                        imageFileName,
                        width,
                        height,
                        froms,
                        0.0d,
                        tos,
                        yMax)

        d.autoSave = false
        d.setMargin(50,50)
        d.color(200,200,200)
        d.drawBorder()
        d.drawRegions()

        d.drawYAxis([0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0], true)

        double minX = targets[0].from
        double maxX = targets[-1].to

        // Line at 1.0 in black
        d.line(displayRegion.from, 1.0, displayRegion.to, 1.0)

        for(target in targets) {
            Region targetRegion = new Region(cnv.chr, target)
            Matrix coverage = getNormalisedCoverage(targetRegion, cnv.sample)
            //                double [] coverageValues = coverage[][0] as double[]


            double [] targetRange = target as double[]
            double [] otherCoverage = coverage.others as double[]
            
            // Can't plot a curve through too few points (exception from loess interpolator)
            if(targetRange.length < 8)
                continue

            LoessInterpolator interp = new LoessInterpolator()
            PolynomialSplineFunction otherMeansFn = interp.interpolate(targetRange, otherCoverage)

            d.color(220,220,220)
            def sds = coverage.collect { c -> c[0] > 1.0d ? 1.0d + c[2] : 1.0d - c[2]}
            d.loess(target, sds)

            d.color("red")
            d.loess(target, coverage.sample, color : { x, y ->
                def otherMean = otherMeansFn.value(x)
                if(otherMean > 30d)
                    "red"
                else
                if(otherMean > 15d)
                    [100,0,100]
                else
                    "blue"
            })
        }

        // Find the genes overlapping the CNV
        d.color("black")
        for(gene in genes) {
            def geneExons = geneDb.getExons(gene).reduce().grep { it.chr == cnv.chr }

            int geneFrom = Math.max(geneExons*.from.min(), displayRegion.from)
            int geneTo = Math.min(geneExons*.to.max(), displayRegion.to)

            d.bar(geneFrom, geneTo, -0.1, gene)
            d.barHeightUp = 0
            d.barHeightDown = 3
            geneExons.eachWithIndex { exon, exonCount ->

                List<Range> intersects = targetRegions.intersect(exon)
                if(!intersects) {
//                    println "Skip exon ${exonCount+1} because no target region overlaps it"
                    return
                }

                //                    Region clipped = exon.intersect(displayRegion)
                //                    println "$exon intersect $displayRegion = $clipped"
                //d.text(clipped.to + (clipped.to - clipped.from)/2, 0.05, "${exonCount+1}")

                for(exonIntersect in intersects) {
                    //                        println "Drawing exon ${exonCount+1}:  $exonIntersect.from-$exonIntersect.to"
                    d.bar(exonIntersect.from, exonIntersect.to, -0.1, null, "${exonCount+1}")
                }
            }
            d.barHeightUp = 5
            d.barHeightDown = 5
        }
        

        double offset = 0d
        for(caller in callers) {

            def callerCnvs = cnvCalls[caller].getOverlaps(cnv).grep { it.extra.sample == cnv.sample }
            if(!callerCnvs) {
                println "No overlapping calls for $caller"
                continue
            }

            d.color(colors[caller])

            println "Caller $caller cnvs are: " + callerCnvs.collect { it.from + " - " + it.to + " (qual=" + it.extra.quality + ")" }
            
            def labels = callerCnvs.collect { (it.extra.quality != null) ? "$caller [" + numberFormat.format(it.extra.quality.toFloat()) + "]" : caller  }
            
            d.bars(callerCnvs*.from, callerCnvs*.to, [1.8 - offset]*callerCnvs.size(), labels)
            offset += 0.08
        }
        
        if(sampleInfo[cnv.sample].files.vcf) {
            if(!vcfs[cnv.sample]) {
                vcfs[cnv.sample] = VCF.parse(sampleInfo[cnv.sample].files.vcf[0]).toBED()
            }
            Regions vcf = vcfs[cnv.sample]
            def variants = null
            synchronized(vcf) {
                variants = targets.collect { vcf.getOverlaps(cnv) }.sum()
            }
            println "Found ${variants.size()} variants for CNV $cnv"
            for(variant in variants*.extra) {
                
                float variantHeight = 0.1f
                d.color("red")
                d.line(variant.pos, yMax, variant.pos, yMax - variantHeight) 
                
                if(variant.sampleDosage(cnv.sample)==1) {
                    int sampleIndex = variant.header.samples.indexOf(cnv.sample)
                    int refReads = variant.getAlleleDepths(0)[sampleIndex]
                    int altReads = variant.getAlleleDepths(1)[sampleIndex]
                    d.color("green")
                    try {
                        d.line(variant.pos, yMax, variant.pos, yMax - variantHeight * ((float)refReads / (refReads+altReads)))
                    }
                    catch(ArithmeticException e) {
                        println "WARNING: Unable to draw variant $variant (refReads=$refReads, altReads=$altReads)"
                    }
                }
            }
        }
        d.save()
        
        if(this.amplicons) {
            d.fileName = d.fileName.replaceAll('.png$','.ac.png')
            drawAmplicons(d, cnv)
            d.save()
            println "Saved " + new File(d.fileName).absolutePath
        }
    }
    
    /**
     * Identify amplicons over the CNV region and draw each amplicon as a bar
     * 
     * @param d
     * @param cnv
     */
    void drawAmplicons(Drawing d, Region cnv) {
        
        d.color("black")
        
        // Find the amplicons overlapping the CNV
        this.ampliconCounts.amplicon = this.amplicons
        
        List<Integer> indices = this.amplicons.findIndexValues { it.overlaps(cnv) }
        
        int sampleIndex = this.ampliconCounts.@names.indexOf(cnv.sample)
        
        Matrix overlapping = this.ampliconCounts[indices]
        
        List<Region> cnvAmplicons = this.amplicons[indices]
        
        println "Found $overlapping.rowDimension amplicons overlapping CNV $cnv for sample $cnv.sample (index=$sampleIndex)"
        
        overlapping.stats = overlapping.collect { Stats.from(it) { value, index ->
            return index != sampleIndex
        }}
        
        println "Means are " + overlapping.stats*.mean
        
        overlapping.relMean = overlapping.collect { values -> values[sampleIndex] / stats.mean }
        
        println "RelMeans are " + overlapping.relMean
        
        if(noAmpliconLabels) {
            d.bars(cnvAmplicons*.from, cnvAmplicons*.to, overlapping.relMean)            
        }
        else {
            d.fontSize(7)
            d.bars(cnvAmplicons*.from, cnvAmplicons*.to, overlapping.relMean, overlapping.collect { values ->
                numberFormat.format(values[sampleIndex]) + " / " + numberFormat.format(stats.mean) + "±" + numberFormat.format(stats.standardDeviation)
            })
            d.fontSize(11)
        }
    }
    
    Matrix getNormalisedCoverage(Region region, String sample) {
        
        Matrix normCovs = normaliseCoverage(region)
        
        int sampleIndex = allSamples.indexOf(sample)
        
        // Get sample relative to mean of other samples
        Matrix rel = normCovs.transformRows { row ->
            
            // Note: we want the coverage for all samples EXCEPT the one that the CNV was called in
            Stats stats = Stats.from(row) { x, index -> index != sampleIndex }
            double effectiveMean = stats.mean > 0.0d?stats.mean : 1
            [row[sampleIndex] / effectiveMean, effectiveMean, stats.standardDeviation / effectiveMean]
        }
        
        rel.@names = ["sample", "others", "sd"]
        
        return rel
    }

    def Matrix normaliseCoverage(Region region) {
        
        double [] sampleMeans = allSamples.collect { means[it]?:1 } as double[]
        
        Matrix sampleCovs = getCoverageMatrix(region)

        double meanOfMeans = Stats.mean(sampleMeans)

        Matrix normCovs = sampleCovs.transform { x, row, col ->
            (x / sampleMeans[col]) * meanOfMeans
        }
        return normCovs
    }
    
    ThreadLocal<Map<String, SAMFileReader>> readers = new ThreadLocal()

    @CompileStatic
    Matrix getCoverageMatrix(Region region) {
        
        Matrix sampleCovs = new Matrix(allSamples.collectEntries { s ->
            //            println "Querying coverage for sample $s"
            
            if(readers.get() == null)
                readers.set([:])
                
            SAMFileReader reader = readers.get().get(s)
            if(reader == null) {
                reader = bams[s].newReader()
                readers.get().put(s, reader)
            }
            
            PileupIterator pileup
            try {
                pileup = bams[s].pileup(reader, region.chr, region.from, region.to)
                def sampleCov = pileup.collect { PileupIterator.Pileup pi -> pi.alignments.size() }
                [s, sampleCov ]
            }
            finally {
                if(pileup != null) {
                    pileup.readIterator.close()
                }
            }
        })
        return sampleCovs
    }
    
    boolean validate(boolean ignoreMissing=false) {
        def missingCovs = this.samples.grep { !sampleInfo[it]?.files?.bam }.unique()
        if(missingCovs) {
            if(ignoreMissing) {
                this.samples.removeAll(missingCovs)
            }
            else {
                err.println "The following samples do not have coverage information provided: $missingCovs"
                return false
            }
        }
        
        return true
    }
    
    /**
     * Parse an option specifying a CNV caller
     * 
     * @param caller    the caller to parse the config for
     * @param opt       a list of the parameters supplied for the given caller
     * @param factory   a closure that creates a parser for the caller
     */
    static void parseCallerOpt(String caller, List<String> opt, Closure factory, Map<String,RangedData> results) {
        for(String cfg in opt) { 
            List<String> parts = cfg.tokenize(":")
            if(parts.size()==1) {
                parts.add(0,caller)
            }
            results[parts[0]] = factory(parts[1]).load()
        }        
    }
    
    
    public static void main(String [] args) {
        
        println "=" * 100
        println "SeeCNV Diagram"
        println "=" * 100
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage: "CNVDiagram <options>")
        cli.with {
            cnvs 'Consolidated CNV summary report', args:1
            region 'A single region to plot, alternative to providing -cnvs', args:1
            xhmm 'CNV calls from XHMM', args:Cli.UNLIMITED
            cnmops 'CNV calls from CN.mops', args:Cli.UNLIMITED
            cfr 'CNV calls from Conifer', args:Cli.UNLIMITED
            angel 'Deletion calls from Angel', args:Cli.UNLIMITED
            ed 'CNV calls from Exome Depth', args:Cli.UNLIMITED
            generic  'CNV calls in BED format, sample in id column', args:Cli.UNLIMITED
            vcf 'VCF files containing variant calls for samples (optional) for annotation', args:Cli.UNLIMITED
            bam 'BAM file, one for each sample', args:Cli.UNLIMITED
            bamDir 'Directory to search for BAM files', args:Cli.UNLIMITED
            targets 'BED file containing target regions', args:1, required:true
            amplicons 'BED file containing amplicons (for HaloPlex data)', args:1
            ampliconcounts 'BED file containing read counts for amplicons (for HaloPlex data)', args:1
            sample 'Sample to export for (specify multiple times or omit for all)', args:Cli.UNLIMITED
            refseq 'RefSeq Genes file from UCSC for annotating genes', args: 1
            t 'Number of threads to use (1)', args:1
            o 'Output file (png format)', args:1, required:true
            w 'Width of output file (pixels)', args:1
            h 'Height of output file (pixels)', args:1
            sampleinfo 'A sample meta data file containing VCF and BAM files (replaces -vcf and -bam)', args:1
            xmean 'Use extended, but more accurate / slower estimator of sample mean read depth'
            litemean 'Use extremely quick, less accurate mean coverage estimator'
            gatkcov 'Use precalculated coverage information from gatk, located in given directory', args:1
            ignoremissing 'Ignore samples that are missing coverage information'
            noamplabels 'Omit drawing read counts over each amplicon'
            chr 'Limit drawing to CNVs overlapping given chromosome (eg: chrX)', args:1
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        Map<String,RangedData> cnvCalls = [:]
           
        if(opts.eds) 
            parseCallerOpt("ed", opts.eds, { new ExomeDepthResults(it) }, cnvCalls)
            
        if(opts.xhmms) 
            parseCallerOpt("xhmm", opts.xhmms, { new XHMMResults(it) }, cnvCalls)
        
        if(opts.cnmopss)
            parseCallerOpt("cnmops", opts.cnmopss, { new CNMopsResults(it) }, cnvCalls)
            
        if(opts.angelss)
            parseCallerOpt("angel", opts.angels, { new AngelResults(it) }, cnvCalls)
            
        if(opts.cfrs)
            parseCallerOpt("cfr", opts.cfrs, { new ConiferResults(it) }, cnvCalls)
            
        if(opts.generics) {
            opts.generics.each { cnvBedFileAndName ->
                
                if(!cnvBedFileAndName.contains(":"))
                    throw new Exception("Please specify -generic argument in the form <id>:<file>")
                    
                def (name,bedFile) = cnvBedFileAndName.split(":")
                def bed = new BED(bedFile, withExtra:true).load()
                cnvCalls[name] = new Regions()
                bed.each { cnv ->
                    cnv.sample = cnv.extra
                    cnv.extra = cnv
                    cnvCalls[name].addRegion(cnv)
                }
            }
        }
        
        if(opts.angel) {
            cnvCalls.angelhmm = new AngelResults(opts.angel).load()
        }
            
        if(cnvCalls.isEmpty() && !opts.region) {
            System.err.println "\nERROR: Please give at least one of the xhmm, ed, cnmops, angel or generic arguments\n"
            System.exit(1)
        }
        
        log.info "BAM Files are $opts.bams"
        
        if(opts.samples) {
            log.info "Diagrams will be drawn only for $opts.samples"
        }
        
        List bamFiles = opts.bams ? opts.bams : []
        
        if(opts.bamDir) {
            new File(opts.bamDir).eachFile { f ->
                if(f.name.endsWith(".bam"))
                    bamFiles << f.absolutePath
            }
        }
        
        Regions targetRegions = new BED(opts.targets).load().reduce()
        
        def sampleInfo = null
        if(opts.sampleinfo)
            sampleInfo = SampleInfo.parse_sample_info(opts.sampleinfo)
        else
            sampleInfo = SampleInfo.fromFiles(bamFiles + (opts.vcfs?:[]))
            
        if(!opts.cnvs && !opts.region) {
            System.err.println "ERROR: please provide either -cnvs or -region"
            System.exit(1)
        }
        
        def cnvs = loadCNVs(opts)
            
        int width = opts.w ? opts.w.toInteger() : 1024
        int height = opts.h ? opts.h.toInteger() : 480
        
        def diagram = new CNVDiagram(cnvs, sampleInfo, cnvCalls, targetRegions, opts.samples?:null)
        
        if(opts.t)
            diagram.concurrency = opts.t.toInteger()
            
        if(opts.gatkcov) {
            diagram.gatkMeanEstimator = opts.gatkcov
            log.info "Using coverage estimates from GATK DepthOfCoverage output"
        }
        else
        if(opts.xmean)
            diagram.extendedMeanEstimator = true
        else
        if(opts.litemean)
            diagram.liteMeanEstimator = true
            
        if(!diagram.validate(opts.ignoremissing)) {
            diagram.err.println "ERROR: one or more arguments was missing or invalid"
            System.exit(1)
        }
        
        if(opts.refseq) {
            diagram.geneDb = new RefGenes(opts.refseq, stripChr: cnvs.numberOfRanges>0 && !cnvs[0].chr.startsWith('chr'))
        }
        
        if(opts.amplicons) {
            diagram.amplicons = new BED(opts.amplicons,withExtra:true).load()
            println "Loaded ${diagram.amplicons.numberOfRanges} amplicons"
        }
       
        if(opts.ampliconcounts) {
            diagram.ampliconCounts = Matrix.load(opts.ampliconcounts)
            println "Loaded read counts for ${diagram.ampliconCounts.rowDimension} amplicons x $diagram.ampliconCounts.columnDimension samples"
        }
        
        diagram.noAmpliconLabels = opts.noamplabels
        
            
        diagram.draw(opts.o, width, height)
    } 

    private static Regions loadCNVs(OptionAccessor opts) {
        Regions cnvs
        if(opts.region) {
            def parts = (opts.region =~ /(.*?):(.*$)/)[0]
            if(!parts)
                err "Please specify the region to draw in the form <sample>:chrX:00000-11111"

            Region cnv =  new Region(parts[2])
            cnv.start = cnv.from
            cnv.end = cnv.to
            cnv.sample = parts[1]
            cnvs = new Regions()
            cnvs.addRegion(cnv)
        }
        else
        if(opts.cnvs.endsWith(".bed")) {
            cnvs = new Regions()
            new BED(opts.cnvs, withExtra:true).load().each { Region r ->
                r.sample = r.extra
                r.extra = r
                cnvs.addRegion(r)
            }
        }
        else
            cnvs = new RangedData(opts.cnvs).load()
            
        if(opts.chr) {
            cnvs = cnvs.grep { it.chr == opts.chr } as Regions
        }
        
        return cnvs
    }
    
    private Regions findLiteMeanRegions() {
        
       int totalRanges = targetRegions.numberOfRanges
       int sampleSize = Math.min(totalRanges, 100)
       Regions result = new Regions()
       Random random = new Random(12345) // fixed seed for reproducibility
       for(int i=0; i<sampleSize;++i) {
           int regionIndex = random.nextInt(totalRanges-1)
           result.addRegion(targetRegions[i])
       }
       
       return result
    }
}
