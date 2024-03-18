import java.awt.BasicStroke;

import java.text.NumberFormat
import java.util.logging.Level

import jsr166y.ForkJoinPool;
import ximmer.CountReadsMeanEstimator
import ximmer.GATKMeanEstimator
import ximmer.MeanEstimator
import ximmer.MultiCovMeanEstimator
import ximmer.PileupMeanEstimator
import ximmer.results.*

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction
import org.codehaus.groovy.runtime.StackTraceUtils

import gngs.*
import gngs.sample.SampleInfo
import gngs.tools.ShearingKmerCounter
import graxxia.Matrix
import graxxia.Stats
import groovy.json.JsonOutput
import groovy.transform.CompileStatic;
import groovy.transform.Memoized
import groovy.util.logging.Log;
import groovy.xml.Namespace;
import groovyx.gpars.GParsPool;
import htsjdk.samtools.SamReader

/**
 * Draws a diagram of the region of each CNV in a provided CNV summary report
 * 
 * @author simon
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
    
    /**
     * Global mean coverage for each sample
     */
    Map<String,Double> means = null
    
    RefGenes geneDb = null
    
    Regions targetRegions = null
    
    FASTA reference = null
    
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
    
    String covjsMeanEstimator = null
    
    boolean redraw = false
    
    double yMax = 2.0d
    
    /**
     * What kind(s) of output to produce
     */
    List<String> writeTypes = ['png']
    
    /**
     * Trying to display too many target regions in plots makes them unreadable.
     * This count limits the number of flanking regions that will be drawn next
     * to each deletion in cases where the affected gene is very large or spans a
     * huge number of target regions
     */
    int maxTargetCount = 20
    
    CNVDiagram(Regions cnvs, Map<String,SampleInfo> sampleInfo, Map<String, RangedData> cnvCalls, Regions targetRegions, List<String> vcfFiles, List<String> samples=null) {
        
        this.cnvs = cnvs;
        
        this.cnvs.each { cnv ->
            if(cnv.sample == null)
                cnv.sample = cnv.extra
        }
        
        this.sampleInfo = sampleInfo
        if(samples) {
            this.samples = samples
        }
        else  {
            this.samples = this.cnvs*.sample.unique()
            log.info "Extracted samples from CNVs: ${this.samples.join(',')}"
        }
        
        allSamples = (sampleInfo.keySet() as List).grep { 
            sampleInfo[it].files.bam
        }
        this.cnvCalls = cnvCalls
        this.targetRegions = targetRegions 

        log.info "There are ${this.cnvs.numberOfRanges} regions, starting with ${this.cnvs[0]}"

        
        if(!targetRegions[0].chr.startsWith('chr')) {
            log.info "Replacing chr on CNVs due to target regions using non-chr prefixed regions"
            this.cnvs = this.cnvs.collect { it.setChr(it.chr.replaceFirst('chr',''));  it } as Regions
        }

        log.info "There are ${this.cnvs.numberOfRanges} regions, starting with ${this.cnvs[0]}"

        
        log.info "Parsing VCFs"
        List<VCF> vcfs = vcfFiles.collect { VCF.parse(it) { cnvs.overlaps(it) } }
        
        Map<VCF, Regions> vcfRegions = vcfs.collectEntries { [ it,  it.toRegions() ] }
        
        this.vcfs = allSamples.collectEntries { String sample ->
            VCF sampleVCF = vcfs.find { sample in it.samples }
            [sample, sampleVCF? vcfRegions[sampleVCF] : null ]
        }
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
            MeanEstimator meanEstimator = createMeanEstimator()
            this.means = meanEstimator.calculateMeans(allSamples)
        }
        log.info "Sample means are $means"
    }

    private MeanEstimator createMeanEstimator() {
        MeanEstimator meanEstimator
        if(gatkMeanEstimator) {
            meanEstimator = new GATKMeanEstimator(this.gatkMeanEstimator)
            log.info "Discovered ${meanEstimator.intervalFiles.size()} samples with GATK coverage information"
        }
        else
        if(covjsMeanEstimator) {
            meanEstimator = new MultiCovMeanEstimator(Utils.reader(covjsMeanEstimator))
        }
        else
        if(extendedMeanEstimator || liteMeanEstimator) {
            
            Regions meanRegions = targetRegions
            if(liteMeanEstimator) {
                // Search for the smallest chromosome that we can that is in the target regions
                meanRegions = findLiteMeanRegions()
            }

            if(extendedMeanEstimator)
                log.info "Using extended mean coverage estimation ..."
            else
                log.info "Using lite mean coverage estimation ..."

            meanEstimator = new PileupMeanEstimator(meanRegions, bams)
        }
        else {
            meanEstimator = new CountReadsMeanEstimator(targetRegions, bams)
        }
        return meanEstimator
    }
    
    void draw(String outputFileBase, int width, int height) {
        
        if(!cnvs.count { 1 }) { // if no CNV calls, short circuit loading all the coverage info
            log.info "No CNVs to draw!"
            return
        }
        
        if(!bams)
            loadBAMs()
        
        if(!means) 
            loadMeans()
        
        def callers = cnvCalls.keySet() as List
        List<String> palette = ["red","green","orange","blue","gray","magenta","yellow","cyan","black"]*10
        def colors = [ callers, palette[0..<callers.size()] ].transpose().collectEntries()
        
        Closure processCnv = { cnv ->
            Utils.time("rendering of CNV $cnv.chr:$cnv.from-$cnv.to", log: log, suppressStartMessage:true ) {
                try {
                    drawCNV(cnv, outputFileBase, width, height, callers, colors)
                }
                catch(Exception e) {
                    StackTraceUtils.sanitize(e)
                    log.log(Level.SEVERE, "Failed to draw CNV $cnv", e)
                }
            }
        }

        Utils.time("Draw all CNVs", log:log) {
            if(concurrency > 1) {
                GParsPool.withPool(concurrency) {
                    cnvs.eachParallel(processCnv)
                }
            }
            else {
                cnvs.each(processCnv)
            }
        }
        
        
        if(coverageCachedReads>0) 
            log.info "Finished drawing.  Executed ${coverageReads}/${coverageCachedReads} coverage read operations (caching hit rate = ${Utils.perc(coverageReads/coverageCachedReads)})"
        else
            log.info "Finished drawing (no coverage reads)"
    }
    
    @CompileStatic
    double round2Digits(double x) {
        Math.round(x*100.0f)/100.0
    }
    
    @CompileStatic
    String normChr(String chr) {
       chr.startsWith('chr') ? chr : 'chr' + chr 
    }

    void drawCNV(Region cnv, String outputFileBase, int width, int height, List callers, Map colors) {
        
        log.info "Drawing cnv $cnv in sample $cnv.sample"
        
        if(this.samples && !(cnv.sample in this.samples)) {
            log.info "Skip CNV $cnv.chr:$cnv.from-$cnv.to for sample $cnv.sample"
            return
        }
        
        log.info "Call $cnv.chr:$cnv.from-$cnv.to for sample $cnv.sample"

        String imageFileName = outputFileBase.replaceAll('.png$','') + "_${normChr(cnv.chr)}_${cnv.from}_${cnv.to}_${cnv.sample}.png"
        String jsonFileName = outputFileBase.replaceAll('.png$','') + "_${normChr(cnv.chr)}_${cnv.from}_${cnv.to}_${cnv.sample}.js"
        
        File jsonFile
        Writer json = new StringWriter()
        if('json' in writeTypes) {
            log.info "Enabling JSON output"
            jsonFile = new File(jsonFileName)
        }
        
        File imageFile = new File(imageFileName)
        if(imageFile.exists() && (imageFile.length() > 0) && (jsonFile != null && jsonFile.exists()) && !redraw) {
            log.info "Skip $imageFileName because file jsonFile.absolutePath already exists"
            return
        }
        
        SAM bam = bams[cnv.sample]

        List<String> genes = geneDb.getGenes(cnv).grep { gene -> geneDb.getExons(gene).overlaps(this.targetRegions) }

        List<IntRange> targets = determineCnvTargets(cnv, genes)
       
        def froms = targets*.from
        def tos = targets*.to
        
        if(froms.size()==0 || tos.size()==0) {
            log.warning  "ERROR: no targets overlapped by $cnv.chr:$cnv.from-$cnv.to"
            return
        }
        
        if(jsonFile)
            json = jsonFile.newWriter() 
        
        json.println("""cnv = //NOJSON\n{ 
     "start" : ${cnv.from},
     "end" : ${cnv.to},
     "targets" : [""")
        
        
        Region displayRegion = new Region(cnv.chr, froms[0]..tos[-1])
        
        boolean drawPNG = 'png' in writeTypes 

        graxxia.Drawing d
        if(drawPNG) 
            d = createDiagramFrame(imageFileName, displayRegion, targets, width, height)

        double minX = targets[0].from
        double maxX = targets[-1].to

        boolean first = true
        for(target in targets) {
            
            if(!first)
                json.println(",")
                
            first = false
            
            outputTargetCoverage(cnv, target, json, d)
        }
        
        json.println("\n    ],")
        json.println("""  "genes" : [""")

        if(d)
            d.color("black")
            
        // For each genes overlapping the CNV
        first = true
        for(gene in genes) {
            
            if(!first)
                json.println(',')
            
            first = false
            
            def geneExons = geneDb.getExons(gene).reduce().grep { it.chr == cnv.chr }

            int geneFrom = displayRegion.from
            if(geneExons)
                geneFrom = Math.max(geneExons*.from.min(), displayRegion.from)

            int geneTo = displayRegion.to
            if(geneExons)
                geneTo = Math.min(geneExons*.to.max(), displayRegion.to)
            
            Map jsonGene = [
                gene: gene,
                start: geneFrom,
                end: geneTo,
                exons: []
            ]
            
            if(d) {
                d.bar(geneFrom, geneTo, -0.1, gene)
                d.barHeightUp = 0
                d.barHeightDown = 3
            }
            
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
                    if(d)
                        d.bar(exonIntersect.from, exonIntersect.to, -0.1, null, "${exonCount+1}")
                    jsonGene.exons << [ 
                        from: exonIntersect.from, 
                        to: exonIntersect.to,
                        number: exonCount+1
                    ]
                }
            }
            
            if(d) {
                d.barHeightUp = 5
                d.barHeightDown = 5
            }
            
            json.print("      " + JsonOutput.toJson(jsonGene))
        }
        
        json.println("\n    ],")
        json.println("""   "callers" : [""")

        double offset = 0d
        first = true
        for(caller in callers) {
            
            def callerCnvs = cnvCalls[caller].getOverlaps(cnv).grep { it.extra.sample == cnv.sample }
            if(!callerCnvs) {
                log.info "No overlapping calls of $cnv for $caller"
                continue
            }
            
            if(!first)
                json.print(',')
            
            first = false

            if(d) 
                d.color(colors[caller])

            log.info "Caller $caller cnvs are: " + callerCnvs.collect { it.from + " - " + it.to + " (qual=" + it.extra.quality + ")" }
            
            List labels = callerCnvs.collect { (it.extra.quality != null) ? "$caller [" + numberFormat.format(it.extra.quality.toFloat()) + "]" : caller  }
            
            Map callerJson = [
                caller: caller,
                calls : callerCnvs.collect {[
                    start: it.from,
                    end: it.to,
                    quality: it.extra.quality?.toFloat()
                ]}
            ]
            
            if(d)
                d.bars(callerCnvs*.from, callerCnvs*.to, [1.8 - offset]*callerCnvs.size(), labels)
            offset += 0.08
            json.print("      " +  JsonOutput.toJson(callerJson))
        }
        
        json.println("\n    ],")
        json.println("""   "variants" : [""")
        if(vcfs[cnv.sample]) {
            drawVariants(json, d, cnv, targets)
        }
        else {
            log.info "Sample $cnv.sample does not have an associated VCF"
        }
        
        if(drawPNG)
            d.save()
        
        if(this.amplicons && d) {
            d.fileName = d.fileName.replaceAll('.png$','.ac.png')
            drawAmplicons(d, cnv)
            d.save()
            println "Saved " + new File(d.fileName).absolutePath
        }
        
        if(json) {
            json.println("    ]\n}")
            json.close()
        }
    }

    void outputTargetCoverage(Region cnv, IntRange target, Writer json, graxxia.Drawing d) {
        
        boolean drawPNG = 'png' in writeTypes 
        
        Region targetRegion = new Region(cnv.chr, target)
        
        // Exclude all other samples that have a CNV call in this region
        // UNLESS that leaves less than 3 samples?
        List<Region> otherCnvs = cnvs.getOverlaps(cnv)*.extra.grep { it.sample != cnv.sample }
        
        List<String> excludeSamples = otherCnvs*.sample.unique()
        
        Matrix coverage = getStandardisedCoverage(targetRegion, cnv.sample, excludeSamples)

        double [] targetRange = target as double[]

        outputTargetJSONCoverage(json, targetRegion, coverage)

        // Can't plot a curve through too few points (exception from loess interpolator)
        if(targetRange.length < 8)
            return 

        if(drawPNG) {
            drawInterpolatedCoverage(d, coverage, target, targetRange)
        }
    }
    
    void drawInterpolatedCoverage(graxxia.Drawing d, Matrix coverage, IntRange target, double [] targetRange) {
        double [] otherCoverage = coverage.others as double[]
        def sds = coverage.collect { c -> c[0] > 1.0d ? 1.0d + c[2] : 1.0d - c[2]}
        LoessInterpolator interp = new LoessInterpolator()
        PolynomialSplineFunction otherMeansFn = interp.interpolate(targetRange, otherCoverage)
        d.color(220,220,220)
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
    
    void drawVariants(Writer json, graxxia.Drawing d, Region cnv, List<IntRange> targets) {
        
        Regions vcf = vcfs[cnv.sample]
            
        def variants = null
        
        synchronized(vcf) {
            variants = targets.collect { vcf.getOverlaps(cnv)*.extra.grep { it.sampleDosage(cnv.sample) > 0 } }.sum()
        }
        
        log.info "Found ${variants.size()} variants for CNV $cnv"
        for(Variant variant in variants) {
            float variantHeight = 0.1f
                
            if(d) {
                d.color("red")
                d.line(variant.pos, yMax, variant.pos, yMax - variantHeight) 
            }
                
            float altFrac = 1.0f
            int dosage = variant.sampleDosage(cnv.sample)
            int sampleIndex = variant.header.samples.indexOf(cnv.sample)
            int altReads = variant.getAlleleDepths(1)[sampleIndex]
            int refReads = variant.getAlleleDepths(0)[sampleIndex]
            if(dosage==1) {
                if(d)
                    d.color("green")
                try {
                    altFrac = ((float)refReads / (refReads+altReads))
                        
                    if(d) 
                        d.line(variant.pos, yMax, variant.pos, yMax - variantHeight * altFrac)
                }
                catch(ArithmeticException e) {
                    log.warning "WARNING: Unable to draw variant $variant (refReads=$refReads, altReads=$altReads)"
                }
            }
                
            Map variantJson = [
                pos: variant.pos,
                dosage: dosage,
                frac: altFrac.isNaN() ? 0.0f : altFrac,
                ref: refReads,
                alt: altReads
            ]
            json.println("      " +  JsonOutput.toJson(variantJson) + ",")
        }
    }
    
    void outputTargetJSONCoverage(Writer json, Region targetRegion, Matrix coverage) {
        Map targetJson = [
            start: targetRegion.from,
            end: targetRegion.to,
            sampleCov: coverage.sample.collect { round2Digits(it) },
            otherCov: coverage.others.collect { Math.round(it) },
            coverageSd: coverage.sd.collect { round2Digits(it) }
        ]
        
        if(this.reference) {
            targetJson.gc = this.reference.gc(targetRegion)
        }
        
        json.print("      " + JsonOutput.toJson(targetJson))
    }
    
    graxxia.Drawing createDiagramFrame(String imageFileName, Region displayRegion, List<IntRange> targets, int width, int height) {
        
        double minX = targets[0].from
        double maxX = targets[-1].to
        
        graxxia.Drawing d = new graxxia.Drawing(
                imageFileName,
                width,
                height,
                targets*.from,
                0.0d,
                targets*.to,
                yMax)

        d.autoSave = false
        d.setMargin(50,50)
        d.color(200,200,200)
        d.drawBorder()
        d.drawRegions()

        d.drawYAxis([0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0], true)
        
        // Line at 1.0 in black
        d.line(displayRegion.from, 1.0, displayRegion.to, 1.0)
        
        return d
    }
    
    /**
     * Decide on a list of target regions to plot for this CNV
     * <p>
     * The list of regions starts as the CNV itself. Then, it is
     * expanded so that any partially overlapped are fully included.
     * Then finally, we want to ensure there is always some context to the
     * left and right of the CNV - so if the CNV is abutting the very
     * edge on either side then the plot regions are expanded one region
     * on each side.
     * 
     * @param cnv   the CNV being plotted (as a Regions object)
     * @param genes a list of genes that should be plotted
     * @return   a list of intervals to be plotted
     */
    List<IntRange> determineCnvTargets(Region rawCNV, List<String> genes) {

        log.info "Genes overlapped by $rawCNV are $genes"

        List targets = null
        
        // We assume later on that the CNV call boundaries line up with the 
        // target region boundaries. To avoid breaking this assumption, 
        // intersect them now
        
        List<Range> cnvTargetIntersect = this.targetRegions.intersect(rawCNV)
        
        Region cnv 
        if(cnvTargetIntersect.isEmpty()) {
            cnv = rawCNV
            log.warning "Error: CNV $rawCNV does not intersect any target regions!"
        }
        else {
            cnv = new Region(rawCNV.chr, cnvTargetIntersect[0].from, cnvTargetIntersect[-1].to)
        }
        
        // If one or more genes is overlapped, display the whole gene(s)
        if(!genes.empty) {
    
            IntRange geneExonsRange = resolveOverlappingGeneExonsRange(genes, cnv)
            
            // Find all the target regions that overlap these exons
            targets = targetRegions.getOverlaps(cnv.chr, geneExonsRange.from, geneExonsRange.to)
            
            log.info "Identified ${targets.size()} overlapping targets for $cnv" 
            
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
        
        if(targets.isEmpty()) {
            log.warning("Failed to identify target region for CNV $cnv (genes=$genes)")
            return []
        }
        
        targets = padTargets(cnv, targets)
        
        return targets
    }
    
    /**
     * @param cnv
     * @param targets
     * @return         the given target ranges padded with 2 extra regions
     *                 if the CNV ends at the target boundary
     */
    List<Range> padTargets(Region cnv, List<Range> targets) {
        
        // If the first region is the same as the start of the CNV, show one region to left as well
        // for context
        int minContext = 50
        if(Math.abs(cnv.from - targets[0].from)<minContext) {
            log.info "Pad targets for $cnv to ensure left context"
            IntRange leftTarget = targetRegions.previousRange(cnv.chr, targets[0].from)
            if(leftTarget) {
                IntRange leftLeftTarget = targetRegions.previousRange(cnv.chr, leftTarget.from)
                if(leftLeftTarget)
                    targets = [leftLeftTarget,leftTarget] + targets
                else
                    targets = [leftTarget] + targets
            }
        }
        
        if(Math.abs(cnv.to - targets[-1].to)<minContext) {
            log.info "Pad targets for $cnv to ensure right context"
            IntRange rightTarget = targetRegions.nextRange(cnv.chr, targets[-1].to)
            if(rightTarget) {
                IntRange rightRightTarget = targetRegions.nextRange(cnv.chr, rightTarget.to)
                if(rightRightTarget)
                    targets = targets + [rightTarget, rightRightTarget]
                else
                    targets = targets + [rightTarget]
            }
        } 
        
        return targets
    }
    
    /**
     * Return a suitable range for plotting that encompasses both exons of the
     * given genes and the CNV specified
     * 
     * @param genes
     * @param cnv
     * @return  range encompassing exons of given genes AND the given CNV
     */
    IntRange resolveOverlappingGeneExonsRange(List<String> genes, Region cnv) {
        
        Regions geneExons = genes.sum { geneDb.getExons(it) }.reduce().grep { it.chr == cnv.chr } as Regions
    
//        log.info "Exons are " + geneExons.collect { it.from + "-" + it.to }.join(",")
            
        int exonsStart = geneExons[0].from 
        int exonsEnd = geneExons[-1].to 
            
        log.info "Overlapping exons of $genes are " + geneExons + " from " + exonsStart + " to " + exonsEnd
            
        // Problem: if CNV starts / ends in target region that is prior to begin / end 
        // of gene then this puts the targets chosen as *smaller* than the CNV
        // itself. So expand to at least include the CNV itself (may be further expanded
        // to include some more padding below)
        if(exonsStart > cnv.from)
            exonsStart = cnv.from
                
        if(exonsEnd < cnv.to)
            exonsEnd = cnv.to
        
        return exonsStart..exonsEnd
    }
    
    /**
     * Identify amplicons over the CNV region and draw each amplicon as a bar
     * 
     * @param d
     * @param cnv
     */
    void drawAmplicons(graxxia.Drawing d, Region cnv) {
        
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
    
    final int minDiagramSamples = 5
    
    /**
     * Compute coverage for the given sample standardised to the distribution
     * of other samples over the specified region.
     * 
     * @param region
     * @param sample
     * @return  Matrix with sample, others and sd columns, and a row for each base
     *          within the given region, containing the standardised coverage
     */
    @CompileStatic
    Matrix getStandardisedCoverage(Region region, String sample, List excludeSamples) {
        
        if(allSamples.size() - excludeSamples.size() < minDiagramSamples)  {
            log.info "Cannot exclude ${excludeSamples.size()}/${allSamples.size()} samples from CNV diagram because it would leave < $minDiagramSamples remaining samples"
            excludeSamples = []
        }
        else {
            if(excludeSamples.size()>0)
                log.info "Excluding ${excludeSamples.size()}/${allSamples.size()} samples from standardised coverage of $sample due to overlapping CNV calls" // + excludeSamples.join(',')
        }
            
        Matrix normCovs = normaliseCoverage(region)
        
        int sampleIndex = allSamples.indexOf(sample)
        
        int [] excludeSampleIndices = (excludeSamples + [sample]).collect { allSamples.indexOf(it) } as int[]
        
        // Get sample relative to mean of other samples
        Matrix rel = normCovs.transformRows { double[] row ->
            // Note: we want the coverage for all samples EXCEPT the one that the CNV was called in
            Stats stats = Stats.from(row) { double x, int index -> !excludeSampleIndices.contains(index) }
            double effectiveMean = stats.mean > 0.0d?stats.mean : 1
            return [row[sampleIndex] / effectiveMean, effectiveMean, stats.standardDeviation / effectiveMean]
        }
        
        rel.@names = ["sample", "others", "sd"]
        
        return rel
    }

    /**
     * Compute a matrix of coverage values relative to each sample's mean
     * for each base of the given region
     * 
     * @param region
     * @return  Matrix with a column for each sample and a row for each base within the specified region
     *          containing the coverage at that position relative to the sample's overall mean coverage
     */
    @CompileStatic
    Matrix normaliseCoverage(Region region) {
        
        double [] sampleWideMeans = allSamples.collect { means[it]?:1 } as double[]
        
        Matrix regionCovs = getCoverageMatrix(region)

        double meanOfMeans = Stats.mean(sampleWideMeans)

        // Divide each value in the region coverage by the corresponding
        // sample's global mean coverage to produce coverage values 
        // that are normalised for the number of reads in the sample
        Matrix normCovs = regionCovs.transform { double x, int row, int col ->
            (x / sampleWideMeans[col]) * meanOfMeans
        }
        return normCovs
    }
    
    ThreadLocal<Map<String, SamReader>> readers = new ThreadLocal()

    /**
     * Returns a matrix containing raw absolute coverage depth for each sample over the 
     * specified region
     * 
     * @param region
     * @return  Matrix with one column per sample, and a row per base position within the
     *          given matrix
     */
    @CompileStatic
    Matrix getCoverageMatrix(Region region) {
        
        Matrix sampleCovs = new Matrix(allSamples.collectEntries { s ->
            ++coverageCachedReads
            short[] sampleCov = readSampleCoverage(s,region.chr, region.from, region.to)
            return [s, sampleCov ]
        })
        return sampleCovs
    }
    
    int coverageCachedReads = 0
    
    int coverageReads = 0
    
    @Memoized(maxCacheSize=200000)
    short[] readSampleCoverage(final String s, String chr, int from, int to) {
        
        double [] kmerFactors = null
        if(this.allKmerFactors)
            kmerFactors = allKmerFactors[allKmerFactors.sample.indexOf(s)]
        
        ++coverageReads

        if(readers.get() == null)
            readers.set([:])
            
        SamReader reader = readers.get().get(s)
        if(reader == null) {
            reader = bams[s].newReader()
            readers.get().put(s, reader)
        }
        
        PileupIterator pileup
        try {
            pileup = bams[s].pileup(reader, chr, from, to)

            short[] sampleCov
            if(kmerFactors == null)  {
                sampleCov = pileup.collect { PileupIterator.Pileup pi -> 
                    pi.alignments.size() 
                } as short[]
            }
            else {
                sampleCov = calculateKmerAdjustedCoverage(kmerFactors, pileup)
            }

            return sampleCov
        }
        finally {
            if(pileup != null) {
                pileup.readIterator.close()
            }
        }
    }
    
    @CompileStatic
    short [] calculateKmerAdjustedCoverage(double [] kmerFactors, PileupIterator pileup) {
        return pileup.collect { PileupIterator.Pileup pi ->
            double cov = 0.0d
            final List<PileupState> alignments = pi.alignments
            final int n = alignments.size()
            for(int i=0; i<n; ++i) {
                int kmerIndex = ShearingKmerCounter.computeKmer(alignments[i].read)
                cov += kmerFactors[kmerIndex]
            }
            return (short)cov
        } as short[]
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
            results[parts[0]] = factory(parts[1])
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
            savvy 'CNV calls from SavvyCNV', args:Cli.UNLIMITED
            cnmops 'CNV calls from CN.mops', args:Cli.UNLIMITED
            cfr 'CNV calls from Conifer', args:Cli.UNLIMITED
            cdx 'CNV calls from CODEX', args:Cli.UNLIMITED
            dfn 'CNV calls from Delfin', args:Cli.UNLIMITED
            angel 'Deletion calls from Angel', args:Cli.UNLIMITED
            ed 'CNV calls from Exome Depth', args:Cli.UNLIMITED
            kmer 'Enables kmer normalisation via supplied kmer profile', args: 1
            generic  'CNV calls in BED format, sample in id column', args:Cli.UNLIMITED
            vcf 'VCF files containing variant calls for samples (optional) for annotation', args:Cli.UNLIMITED
            bam 'BAM file, one for each sample', args:Cli.UNLIMITED
            bamDir 'Directory to search for BAM files', args:Cli.UNLIMITED
            targets 'BED file containing target regions', args:1, required:true
            amplicons 'BED file containing amplicons (for HaloPlex data)', args:1
            ampliconcounts 'BED file containing read counts for amplicons (for HaloPlex data)', args:1
            sample 'Sample to export for (specify multiple times or omit for all)', args:Cli.UNLIMITED
            refseq 'RefSeq Genes file from UCSC for annotating genes', args: 1
            ref 'Reference sequence. If provided, GC content will be annotated', args:1
            t 'Number of threads to use (1)', args:1
            o 'Output file (png format)', args:1, required:true
            w 'Width of output file (pixels)', args:1
            h 'Height of output file (pixels)', args:1
            sampleinfo 'A sample meta data file containing VCF and BAM files (replaces -vcf and -bam)', args:1
            xmean 'Use extended, but more accurate / slower estimator of sample mean read depth'
            litemean 'Use extremely quick, less accurate mean coverage estimator'
            gatkcov 'Use precalculated coverage information from gatk, located in given directory', args:1
            covjs 'Use precalculated mean coverage information from MultiCov js output', args:1
            ignoremissing 'Ignore samples that are missing coverage information'
            noamplabels 'Omit drawing read counts over each amplicon'
            redraw 'Redraw image even if file already exists'
            chr 'Limit drawing to CNVs overlapping given chromosome (eg: chrX)', args:1
            json 'Include JSON in output'
            nopng 'Suppress writing of PNG files'
            autoFilter 'Auto-filter the CNVs to be drawn to include only those likely to be of interest (ie: low population frequency, non-artefacts)'
            autoFilterSampleCount 'Max sample count to apply when auto-filtering CNVs', args:1, type: Integer
            autoFilterMaxFreq 'Max frequency to allow when applying auto-filtering of CNVs', args: 1, type: Double
        }

        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        Regions targetRegions = new BED(opts.targets).load().reduce()

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
            
        if(opts.cdxs)
            parseCallerOpt("cdx", opts.cdxs, { new CodexResults(it) }, cnvCalls)
             
        if(opts.dfns)
            parseCallerOpt("dfn", opts.dfns, { new DelfinResults(it) }, cnvCalls)

        if(opts.savvys)
            parseCallerOpt("savvys", opts.savvys, { new SavvyCNVResults(it, targetRegions) }, cnvCalls)

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
            System.err.println "\nERROR: Please give at least one of the xhmm, ed, cnmops, savvy, angel or generic arguments\n"
            System.exit(1)
        }
        
        log.info "BAM Files are \n\n${opts.bams.join('\n')}"
        
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
        
       
        def sampleInfo = null
        if(opts.sampleinfo) {
            log.info "Parsing sample info from $opts.sampleinfo"
            sampleInfo = SampleInfo.parse_sample_info(opts.sampleinfo)
        }
        else {
            sampleInfo = SampleInfo.fromFiles(bamFiles + (opts.vcfs?:[]))
            log.info "Extracting sample info from bam files: \n${sampleInfo*.key.join('\n')}"
        }
            
        if(!opts.cnvs && !opts.region) {
            System.err.println "ERROR: please provide either -cnvs or -region"
            System.exit(1)
        }
        
        def cnvs = loadCNVs(opts)
            
        int width = opts.w ? opts.w.toInteger() : 1024
        int height = opts.h ? opts.h.toInteger() : 480
        
        def diagram = new CNVDiagram(cnvs, sampleInfo, cnvCalls, targetRegions, opts.vcfs?:[], opts.samples?:null)
        
        if(opts.t)
            diagram.concurrency = opts.t.toInteger()
            
        if(opts.gatkcov) {
            diagram.gatkMeanEstimator = opts.gatkcov
            log.info "Using coverage estimates from GATK DepthOfCoverage output"
        }
        else
        if(opts.covjs) {
            diagram.covjsMeanEstimator = opts.covjs
            log.info "Using coverage estimates from MultiCov output"
        }
        else
        if(opts.xmean)
            diagram.extendedMeanEstimator = true
        else
        if(opts.litemean)
            diagram.liteMeanEstimator = true
            
        if(opts.ref)
            diagram.reference = new FASTA(opts.ref)
            
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
        
        if(opts.redraw)
            diagram.redraw = true
        
        if(opts.json)
            diagram.writeTypes << 'json'
            
        if(opts.nopng)
            diagram.writeTypes.remove('png')
            
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
            String altChr = 'chr' + opts.chr
            cnvs = cnvs.grep { it.chr == opts.chr || it.chr == altChr } as Regions
        }
        
        if(opts.autoFilter) {
            if(opts.region || opts.cnvs.endsWith('.bed'))
                throw new IllegalArgumentException("The -autoFilter option can only be used when CNVs are loaded using the -cnvs option")

            final int maxSampleCount = (opts.autoFilterSampleCount?:3)
            final double maxFreq = (opts.autoFilterMaxFreq?:0.2)
            
            log.info "Filtering CNVs to draw to sample count < $maxSampleCount, population freq < $maxFreq, and dups with > 1 caller"
            
            int oldCNVCount = cnvs.numberOfRanges

            // Filter high freq, high within-batch calls and dups called by a single caller
            cnvs = cnvs.grep { 
                (it.sampleCount < maxSampleCount) &&  // Do not include if found in more than (default 3) samples in the batch
                (it.DDDFreq<maxFreq) && (it.DGVFreq<maxFreq) && (it.GMDFreq<maxFreq) &&  // Do not include if spanning population freq too high 
                (it.type.contains('DEL') || (it.count > 1))  // Only include non-deletions if count > 1
            } as Regions
            
            log.info "Filtering reduced CNVs from $oldCNVCount to ${cnvs.numberOfRanges} (" + Utils.perc(cnvs.numberOfRanges/(oldCNVCount+1)) + '%)'
        }
        
        return cnvs
    }
    
    Matrix allKmerFactors
    
    private void initKmerProfile() {
        
       log.info "Calculating kmer weightings based on kmer profile: $opts.kmer"
       Matrix kmerCounts = Matrix.load(opts.kmer) 
       Matrix norm = kmerCounts.normaliseRows().normaliseColumns()
       this.allKmerFactors = norm.transform { double value ->
           if(value != 0) {
               return 1.0d / value
           }
           else {
               return 1.0d // multiplier will be 1.0d when factor cannot be calculated
           }
       }.fillna(1.0)
       
       log.info "Calculated kmer profile matrix:\n" + this.allKmerFactors
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
