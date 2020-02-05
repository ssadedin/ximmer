// vim: sw=4 expandtab cindent ts=4
import org.codehaus.groovy.runtime.StackTraceUtils;

import gngs.*
import graxxia.Stats
import groovy.json.JsonOutput
import groovy.text.SimpleTemplateEngine
import groovy.transform.CompileStatic
import groovy.util.logging.Log;

import ximmer.results.*
import ximmer.*

/**
 * Reads results from any number of CNV callers and combines them together into
 * a consolidated report in TSV and HTML format.
 * <p> 
 * @author Simon Sadedin
 */
@Log
class SummarizeCNVs {
    
    /**
     * CNV callers to summarize results from
     */
    Map<String,RangedData> results = [:]
    
    /**
     * Targeted regions by assay, with id in column 4 for annotation (typically gene)
     */
    Regions targetRegions
    
    /**
     * Required for annotation of spanning population CNVs
     */
    String dgvFile = null
    
    /**
     * VCFs containing variants for samples for CNV calls (optional)
     */
    List<VCF> vcfs = []
    
    /**
     * Thresholds for filtering CNVs by quality when writing out CNVs
     */
    Map<String, Float> qualityFilters = [:]
    
    /**
     * Locations of variants overlapping CNV calls 
     */
    Map<String, Regions> variants = [:]
    
    List<String> samples = []
    
    TargetedCNVAnnotator cnvAnnotator = null 
    
    RefGenes refGenes = null
    
    String idMask = null
    
    Set<String> filterToGenes = null
    
    Set<String> excludeGenes = null
    
    /**
     * Global gene categories (applied to all samples without a gene list)
     */
    Map<String, Integer> geneCategories = [:]
    
    /**
     * Per gene list gene categories (applied when a sample has been assigned specifically 
     * to a gene list)
     */
    Map<String, Map<String, Integer>> genelistCategories
    
    /**
     * If genelists have been assigned per sample, they are mapped here
     */
    Map<String, String> sampleToGenelist = [:]
    
    Map<String,String> genelists
    
    gngs.VCFIndex gnomad
    
    /**
     * The gene category below which results will be filtered out. If zero is set,
     * genes with no category will be included. Otherwise, a CNV must be both assigned
     * a category and that category must also exceed or equal the minimum category
     * to avoid being filtered out.
     */
    int minimumCategory = 0
    
    /**
     * The fraction of mutual overlap at which calls from differnet CNV callers will be
     * counted as supporting an overall larger call
     */
    double mergeOverlapThreshold = 0.5
    
    /**
     * Method to decide if two overlapping CNVs should be counted as the same CNV
     */
    OverlapCriteria overlapCriteria
    
    /**
     * Regions to exclude CNVs from when they overlap by at least 50%
     */
    Regions excludeRegions50
    
    int countExcludedByRegionOverlap = 0
    
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
            results[parts[0]] = factory(*parts[1..-1])
            //results[parts[0]] = factory(parts[1]).load()
        }        
    }
    
    static void main(String [] args) {
        
        println XimmerBanner.banner("Ximmer CNV Summarizer")
        
        Cli cli = new Cli(usage: "SummarizeCNVs <options>")
        
        cli.with {
            ed 'ExomeDepth results', args:Cli.UNLIMITED
            xhmm 'XHMM results', args:Cli.UNLIMITED
            cnmops 'CN Mops results', args:Cli.UNLIMITED
            cfr 'Conifer results', args:Cli.UNLIMITED
            angel 'Angel results', args:Cli.UNLIMITED
            ex 'Excavator results', args:Cli.UNLIMITED
            cdx 'CODEX results', args:Cli.UNLIMITED
            px 'Parallax results', args:Cli.UNLIMITED
            cnvn 'CNVNator results', args:Cli.UNLIMITED
            canv 'Canvas results', args:Cli.UNLIMITED
            delly 'Delly results', args:Cli.UNLIMITED
            lumpy 'Lumpy results', args:Cli.UNLIMITED
            truth 'Postiive control CNVs', args:1
            vcf 'VCF file containing variants for a sample in results', args:Cli.UNLIMITED
            target 'Target regions with id for each region to annotate', args:1, required:true
            chr 'Process CNVs only for given chromosome', args:1
            report 'Template to use for creating HTML report', args: 1
            dgv 'Path to file containing CNVs from database of genomic variants (DGV)', args: 1
            ddd 'Path to file containing CNVs from Decipher (DDD)', args: 1
            gmd 'Path to gnomAD sites VCF to annotate frequencies from gnomAD', args:1
            samples 'Samples to export', args:1
            samplemap 'Map samples to gene lists in a two column, tab separated file', args:1
            quality 'Filtering by quality score, in form caller:quality...', args:Cli.UNLIMITED
            bam 'BAM file for a sample (if provided, used to customize IGV access)', args:Cli.UNLIMITED
            bampath 'Path to BAM files to enable access using IGV (can be URL)', args:1
            name 'Name for report (displayed in header)', args:1
            refgene 'Path to RefGene file downloaded from UCSC to annotate genes (optional)', args:1
            tsv 'Write consolidated CNV report in tab separated format to <file>', args:1
            json 'Write consolidated CNV report in json separated format to <file>', args:1
            imgpath 'Additional path to add for images', args:1
            genome 'Genome build to use when annotating CNVs', args:1
            idmask 'Mask to apply to sample ids for presentation in report', args:1
            genefilter 'Optional file of genes to filter CNVs to', args:1
            exgenes 'Optional file of genes to exclude from output', args:1
            genelist 'Define a named gene list, to be passed through to report', args:Cli.UNLIMITED
            mincat 'Set a category below which CNVs will be excluded from results' , args: 1
            mergefrac 'The fraction of mutual overlap at which CNVs should be merged to a single call', args:1
            gnomad 'gnomAD VCF to provide population frequency annotations', args:1
            x50 'Exclude CNVs that overlap these regions (bed file) by more than 50%', args:1, required: false
            mergeby 'one of spanbp or sharedtargets', args: 1, required: false
            o 'Output file name', args:1
        }
        
        Utils.configureSimpleLogging()
        
        log.info "Starting ...."
        
        def opts = cli.parse(args)
        if(!opts)
            System.exit(1)
        
        Map<String,RangedData> results = [:]
        
        if(opts.eds) 
            parseCallerOpt("ed", opts.eds, { new ExomeDepthResults(it) }, results)
            
        if(opts.xhmms) 
            parseCallerOpt("xhmm", opts.xhmms, { new XHMMResults(it) }, results)
        
        if(opts.cnmopss)
            parseCallerOpt("cnmops", opts.cnmopss, { new CNMopsResults(it) }, results)
            
        if(opts.angelss)
            parseCallerOpt("angel", opts.angels, { new AngelResults(it) }, results)
            
        if(opts.cfrs)
            parseCallerOpt("cfr", opts.cfrs, { new ConiferResults(it) }, results)
            
		if(opts.cdxs)
            parseCallerOpt("cdx", opts.cdxs, { new CodexResults(it) }, results) 
            
		if(opts.pxs)
            parseCallerOpt("px", opts.pxs, { fileName -> new ParallaxResults(fileName) }, results) 
            
		if(opts.cnvns)
            parseCallerOpt("cnvn", opts.cnvns, { fileName -> new CNVNatorResults(fileName) }, results) 
			
		if(opts.delly)
            parseCallerOpt("delly", opts.dellys, { fileName -> new DellyResults(fileName) }, results) 
            
		if(opts.canv)
            parseCallerOpt("canvas", opts.canvs, { fileName -> new CanvasResults(fileName) }, results) 
            
		if(opts.lumpy)
            parseCallerOpt("lumpy", opts.lumpys, { fileName -> new LumpyResults(fileName) }, results) 
                
        if(opts.truth) 
            results.truth = new PositiveControlResults(opts.truth).load() 
            
        List exportSamples = null
        if(opts.samples)
            exportSamples= opts.samples.split(",")*.trim() as List    

        String outputFileName = opts.report ? new File(opts.report).name : 'cnv_report.html'
        if(opts.o) 
            outputFileName = opts.o
            
        String reportName
        if(opts.name)
            reportName = opts.name

        Map qualityFilters = [:]
        if(opts.qualitys) {
            qualityFilters = opts.qualitys.collectEntries { qs ->
                def parts = qs.split(":")
                if(parts.size()!=2)
                    throw new IllegalArgumentException("Quality filtering string should be in form <caller>:<quality threshold>")
                    
                log.info "Set quality filter for ${parts[0]} to ${parts[1]}"
                [ parts[0], parts[1].toFloat() ]
            }
        }
        
        List<VCF> vcfList = parseVCFs(opts, results)
		
        Regions target = new BED(opts.target, withExtra:true).load()
        
        log.info "Loaded ${Utils.humanBp(target.size())} target regions"
        
        try {
            
            RefGenes refGenes = (!opts.refgene || "download" == opts.refgene) ? RefGenes.download(opts.genome?:"hg19") : new RefGenes(opts.refgene)  
            log.info "Initialised RefGene annotations"
            Set<String> filterToGenes = null
            if(opts.genefilter) {
               filterToGenes = new File(opts.genefilter).readLines()*.trim() as Set
            }
            
            Set<String> excludeGenes = null
            if(opts.exgenes) {
               excludeGenes = new File(opts.exgenes).readLines()*.trim() as Set
            }
            
            SummarizeCNVs summarizer = new SummarizeCNVs(
                results:results, 
                targetRegions:target, 
                vcfs:vcfList, 
                qualityFilters: qualityFilters,
                excludeGenes:excludeGenes,
                filterToGenes: filterToGenes
            )
            
            initFrequencyAnnotator(opts, target, summarizer)
            
            if(opts.refgene == "download") {
                summarizer.refGenes = refGenes
            }
            else
            if(opts.refgene) {
                summarizer.refGenes = new RefGenes(opts.refgene)
            }
            
            if(opts.idmask) {
                summarizer.idMask = opts.idmask
            }
            
            if(opts.genelist) 
                summarizer.genelists = parseGenelistDefinitions(opts)
                
            if(opts.mincat) 
                summarizer.minimumCategory = opts.mincat.toInteger()
            
            if(opts.samplemap) {
                summarizer.sampleToGenelist = new File(opts.samplemap).readLines()*.tokenize('\t').collectEntries()
                log.info "Read ${summarizer.sampleToGenelist} sample-genelist assignments from $opts.samplemap"
            }
            
            if(opts.x50) {
                summarizer.excludeRegions50 = new BED(opts.x50).load()
                log.info "Regions from $opts.x50 will be excluded (${Utils.humanBp(summarizer.excludeRegions50.size())})"
            }
             
            Regions cnvs = summarizer.run(exportSamples)
			
            if(opts.tsv) {
                summarizer.writeTSV(cnvs, opts.tsv)
            }
            
            if(opts.json) {
                summarizer.writeJSON(cnvs, opts.json)
            }
            
            if(opts.mergefrac) {
                summarizer.mergeOverlapThreshold = opts.mergefrac.toDouble()
            }
            
            if(opts.mergeby) {
                initOverlapCriteria(opts, summarizer)
            }
            
            // If there is a truth set available, for each CNV set,
            // set the true CNVs on the set
            if(results.truth) {
                results.each { caller, calls ->
                    if(caller != "truth") {
                        calls.truth = results.truth
                    }
                }
            }
            
            File outputDirFile = new File(outputFileName).absoluteFile.parentFile
            summarizer.writeCallerJSON(results, new File(outputDirFile, 'cnv_calls.js'))
            
            summarizer.log.info "Writing output to " + outputFileName
            summarizer.writeReport(
                cnvs, 
                reportName, 
                outputFileName, 
                opts.report ?: "cnv_report.html", 
                false, 
                opts.bams ?: [], 
                opts.bampath,
                opts.imgpath ?: "")
        }
        catch(Exception e) {
            StackTraceUtils.sanitize(e)
            e.printStackTrace()
            System.exit(1)
        }
    }
    
    static private void initOverlapCriteria(OptionAccessor opts, SummarizeCNVs summarizer) {
        if(opts.mergeby == 'spanbp') {
            summarizer.overlapCriteria = new OverlapBySpan(minimumFraction:summarizer.mergeOverlapThreshold)
        }
        else
        if(opts.mergeby == 'sharedtargets') {
            summarizer.overlapCriteria = 
                new OverlapByFractionOfTargetRegions(targetRegions:summarizer.targetRegions, minimumFraction:summarizer.mergeOverlapThreshold)
        }
        else {
            throw new IllegalArgumentException("The mergeby criteria specified ($opts.mergeby) is not valid")
        }
    }

    private static initFrequencyAnnotator(OptionAccessor opts, Regions target, SummarizeCNVs summarizer) {
        if(opts.dgv || opts.ddd || opts.gmd) {

            Map databases = [:]
            if(opts.ddd) {
                log.info "Loading Decipher annotations from $opts.ddd ..."
                databases["DDD"] = new DecipherCNVs(opts.ddd).parse()
            }

            if(opts.dgv) {
                log.info "Loading DGV annotations from $opts.dgv ..."
                databases["DGV"] = new DGV(opts.dgv).parse()
            }
            
            if(opts.gmd) {
                log.info "Adding gnomAD annotation source from $opts.gmd"
                databases["GMD"] = new GnomADCNVDatabase(new VCFIndex(opts.gmd))
            }

            summarizer.cnvAnnotator = new TargetedCNVAnnotator(databases, target)
        }
    }
    
    static Map<String,File> parseGenelistDefinitions(OptionAccessor opts) {
        def result = opts.genelists.collectEntries { genelistDefinition -> // name:file
                    genelistDefinition.tokenize('=')
        }
        log.info "Gene lists: " + result
        return result
    }
    
    static List<VCF> parseVCFs(OptionAccessor opts, Map<String,RangedData> results) {
        Regions mergedCalls = results*.value.inject(new Regions()) { Regions regions, RangedData calls  ->
            calls.each { regions.addRegion(it) }
            regions
        }.reduce()
        
        List<VCF> vcfList = opts.vcfs ? opts.vcfs.collect { 
            VCF.parse(it) { Variant v ->
                mergedCalls.overlaps(v)
            }
        } : []
        
        return vcfList
    }
    
    /**
     * Write all the CNV calls as JSON
     * 
     * @param results
     * @param outputFile
     */
    void writeCallerJSON(Map<String, CNVResults> results, File outputFile) {
        outputFile.text = 'var cnv_calls = {\n' + results.collect { String caller, CNVResults calls -> 
            /"$caller" : / + calls.toJson(this.cnvAnnotator)
        }.join(',\n') + '\n}\n'
    }
    
    Regions run(List exportSamples) {
        
        if(this.genelists)
            this.parseGenelists()
        
        if(this.dgvFile && !this.cnvAnnotator)
            this.cnvAnnotator = new TargetedCNVAnnotator(targetRegions, dgvFile)
            
        if(this.overlapCriteria == null)
            this.overlapCriteria = new OverlapBySpan(minimumFraction: this.mergeOverlapThreshold)
        
        extractSamples()
        
        // Try to find a VCF file for each sample
        variants = samples.collectEntries { String sample ->
            [sample, vcfs.find { sample in it.samples }?.toRegions()]
        }

        if(exportSamples == null)
            exportSamples = this.samples

        Regions results = new Regions()
        for(s in exportSamples) {
            Regions sampleCnvs = extractMergedSampleCnvs(s)
            for(cnv in sampleCnvs) {
                log.info "Merging CNV $cnv for sample '$s' type = $cnv.type"
                if(excludeRegions50 && excludeRegions50.intersect(cnv)*.size()?.sum()?.div(cnv.size()) > 0.5d) {
                    ++countExcludedByRegionOverlap
                }
                else {
                    results.addRegion(cnv)
                }
            }
        }
        
        if(countExcludedByRegionOverlap)
            log.info "$countExcludedByRegionOverlap cnvs were excluded because they overlap the exclusion bed file"
        
        // Second phase annotation (annotations that depend on the annotations 
        // created in 1st phase
        metaAnnotate(results)
        
        return results
    }
    
    /**
     * Parse files that list gene symbols with optional tab separated priority
     */
    void parseGenelists() {
        
        // Global categories
        this.geneCategories = [:]
        
        // Per gene-list categories
        this.genelistCategories = [:]
        this.genelists.each { name, file -> 
            Map glCats = genelistCategories[name] = [:]
            new File(file).readLines()*.tokenize('\t').collect { [it[0], it[1]?.toInteger()?:1]}.each { gene, cat ->
               geneCategories[gene] = Math.max(cat, geneCategories.get(gene,cat))
               glCats[gene] = cat
            }
            log.info "Read ${glCats.size()} genes from $file"
        }
    }
    
    void metaAnnotate(Regions results) {
        for(Region cnv in results) {
            cnv.targets = this.targetRegions.getOverlaps(cnv).size()
            cnv.samples = results.grep { it.overlaps(cnv) }*.sample.unique()
            cnv.sampleCount = cnv.samples.size()
            cnv.sampleFreq = (cnv.sampleCount / (double)samples.size())
            
//            log.info "Samples for $cnv = $cnv.samples (count = $cnv.sampleCount)"
        }
    }
    
    void writeTSV(Regions cnvs, String fileName) {
        
        log.info "Writing TSV report to $fileName"
        
        List<String> cnvCallers = results.keySet() as List
        Map<String,String> anno_types = [ "DEL" : "LOSS", "DUP" : "GAIN" ]
        
        List annotations = null
        if(cnvAnnotator) {
             annotations = cnvs.collect { cnv -> 
                cnvAnnotator.annotate(new Region(cnv.chr, cnv.from..cnv.to), anno_types[cnv.type]) 
             }
        }
        
        new File(fileName).withWriter { w ->
            
            List<String> dbIds = cnvAnnotator ? cnvAnnotator.cnvDatabases*.key : []
            
            w.println((["chr","start","end","sample","genes", "type","count","stotal","sampleCount","sampleFreq"] + 
                       dbIds.collect{ dbId -> [dbId,dbId+"Freq"]}.sum() +
                       cnvCallers + 
                       cnvCallers.collect { it+"_qual" }).join("\t"))
            
            cnvs.eachWithIndex { cnv, i ->
                
                List frequencyInfo = dbIds.collect { dbId -> CNVFrequency freqInfo = annotations[i][dbId];  [freqInfo.spanning.size(), freqInfo.spanningFreq] }.sum() 
                
                List line = [
                    cnv.chr, 
                    cnv.from,
                    cnv.to, 
                    cnv.sample,
                    cnv.genes,
                    cnv.type, 
                    cnv.count, 
                    cnv.stotal, 
                    cnv.sampleCount,
                    cnv.sampleFreq
                ] + frequencyInfo +
                cnvCallers.collect { caller ->
                    cnv[caller].best ? "TRUE" : "FALSE"
                }  + cnvCallers.collect { caller ->
                    cnv[caller].best ? cnv[caller].best.quality : 0
                } 
                
                w.println line.join("\t")
            }
        }
    }
    
    void writeJSON(Regions cnvs, String fileName) {
        log.info "Writing JSON report to $fileName"
        new File(fileName).withWriter { w ->
            writeJSON(cnvs, w)
        }
    } 
    
    final Map<String,String> anno_types = [ "DEL" : "LOSS", "DUP" : "GAIN" ]
    
    final static List<String> DEFAULT_JS_COLUMNS = 
        ["chr","start","end","targets","sample","genes","category", "type","count","stotal","sampleCount","sampleFreq"]
    
    void writeJSON(Regions cnvs, Writer w) {
        
        List<String> cnvCallers = results.keySet() as List
        
        List<String> dbIds = cnvAnnotator ? cnvAnnotator.cnvDatabases*.key : []
            
        List<String> columnNames = computeColumns(cnvCallers, dbIds)

        w.println('[')
            
        cnvs.eachWithIndex { Region cnv, int i ->
                
            if(i>0)
                w.println(',')
                    
            Map cnvData = cnvToMap(cnvCallers, dbIds, columnNames, cnv)
			
            w.print(JsonOutput.toJson(cnvData))
        }
        w.println('\n]')        
    }
    
    List<String> computeColumns(List<String> cnvCallers, List<String> dbIds) {
       return DEFAULT_JS_COLUMNS + 
           (refGenes?['cds']:[]) +
           (dbIds.collect { String dbId -> [dbId, dbId + 'Freq'] }.sum()?:[]) +
           cnvCallers + cnvCallers.collect { it+"_qual" } + ['calls']
    }
    
    Map<String, Object> cnvToMap(List<String> cnvCallers, List<String> dbIds, List<String> columnNames, Region cnv) {
        
        List cdsInfo = this.refGenes ? [cnv.cdsOverlap] : []
        
        List frequencyInfo = computeCNVFrequencyInfo(cnv, dbIds)
        
        List callerFlags =  cnvCallers.collect { caller ->
            cnv[caller].best ? "TRUE" : "FALSE"
        }
        
        Map calls = [:]
        for(String caller in cnvCallers) {
            if(cnv[caller]) {
                calls[caller] = cnv[caller].all.collect { [it.from, it.to, it.quality] }
            }
        }
        
        List line = [
            cnv.chr, 
            cnv.from,
            cnv.to, 
            cnv.targets,
            cnv.sample,
            cnv.genes,
            cnv.category,
            cnv.type, 
            cnv.count, 
            cnv.stotal, 
            cnv.sampleCount,
            cnv.sampleFreq
        ] + cdsInfo + frequencyInfo +
        cnvCallers.collect { caller ->
            cnv[caller].best ? "TRUE" : "FALSE"
        }  + cnvCallers.collect { caller ->
            cnv[caller].best ? cnv[caller].best.quality : 0
        } + [calls]
               
		Map data = [columnNames,line].transpose().collectEntries()
        
        return data
    }

    private List computeCNVFrequencyInfo(Region cnv, List<String> dbIds) {
        
        String annoType = anno_types[cnv.type]
        
        Region cnvRegion = new Region(cnv.chr, cnv.from..cnv.to)

        Map<String,CNVFrequency> freqInfos = [:]
        if(cnvAnnotator)
            freqInfos = cnvAnnotator.annotate(cnvRegion, annoType)

        List frequencyInfo = dbIds.collect { dbId ->
            CNVFrequency freqInfo = freqInfos[dbId];
            
            assert freqInfo != null : "Annotations from database $dbId were not applied by the configured annotator"

//            if(freqInfo.spanning.size()>0) {
//                log.info "One or more spanning CNVs from $dbId ws identified for $cnv"
//            }

            [freqInfo.spanning.size(), freqInfo.spanningFreq]
        }.sum()?:[]
        return frequencyInfo
    }
    
    @CompileStatic
    void writeReport(Regions cnvs, 
                     String name, 
                     String fileName, 
                     String reportTemplate="cnv_report.html", 
                     boolean inlineJs=true, 
                     List bamFiles = [], 
                     def bamFilePath=false, 
                     String imgpath="") {
        
        log.info "Writing ${cnvs.numberOfRanges} CNVs to report using template: " + reportTemplate
        
        File outputFile = new File(fileName).absoluteFile
        log.info "Output path = " + outputFile.absolutePath
        
        SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
        String jsFileName = new File(reportTemplate).name.replaceAll('\\.html$','\\.js')
        InputStream templateStream 
        String jsCode = null
        
        File cnvReportFile = new File(reportTemplate)
        
        HTMLAssetSource assetSource
        
        // Avoid using the output file as a template if it happens to exist!
        if(cnvReportFile.exists() && (cnvReportFile.canonicalPath != outputFile.canonicalPath)) {
            File templateParentDir = new File(reportTemplate).absoluteFile.parentFile
            assetSource = new HTMLFileAssetSource(templateParentDir)
            templateStream = new File(reportTemplate).newInputStream()
        }
        else {
            assetSource = new HTMLClassloaderAssetSource()
            templateStream = getClass().classLoader.getResourceAsStream(reportTemplate)
        }
        
        File outputDir = outputFile.parentFile
        HTMLAssets assets = 
            new HTMLAssets(assetSource, outputDir) \
                  << new HTMLAsset(
                    source: 'cnv_report.js',
                    name: 'cnv_report.js'
                ) << new HTMLAsset(
                    source: 'cnv_diagram.js'
                ) << new HTMLAsset(
                    source: 'cnv_report.css'
                ) << new HTMLAsset(
                    source: 'vue.js'
                ) << new HTMLAsset(
                    source: 'require.js'
                ) 
                
                ;
                                
        String renderedAssets = assets.render()
        
        templateStream.withReader { r ->
            new File(fileName).withWriter { w ->
                templateEngine.createTemplate(r).make(
                    imgpath: imgpath,
                    cnvs : cnvs,
                    batch_name : name,
                    cnv_callers : results.keySet() as List,
                    types : ['DUP','DEL'],
                    reportSamples : false,
                    cnvAnnotator : cnvAnnotator,
                    js : renderedAssets,
                    bam_files : bamFiles,
                    bam_file_path : bamFilePath,
                    idMask: idMask,
                    geneCategories: geneCategories
                ).writeTo(w)
            }
        }
    }
    
    void extractSamples() {
        this.samples = results.collect { e ->
            e.value*.sample
        }.flatten().unique()
        
        log.info "Extracted samples from calls: $samples"
    }
    
    Regions extractMergedSampleCnvs(String sample) {
        
        // First accumulate CNVs from all the callers into one merged object
        Regions sampleCNVs = extractSampleIndividualCnvs(sample)
        
//        Regions flattened = merged.reduce()
       
        CNVMerger cnvMerger = new CNVMerger(sampleCNVs, this.overlapCriteria)
        Regions flattened = cnvMerger.merge()
        
        log.info "Sample $sample has ${flattened.numberOfRanges} CNVs"
        
        List callers = results.keySet() as List

        Regions result = new Regions()

        // Now annotate each flattened CNV with the callers that found it, and the 
        // confidence of the highest confidence call from each caller
        for(Region cnv in flattened) {
            
            annotateCNV(sample, callers, cnv)
            
            if(!isFiltered(callers, cnv))
                result.addRegion(cnv)
            log.info "Annotated CNV $cnv (${cnv.hashCode()}) for $sample of type $cnv.type"
        }
        
        for(Region cnv in result) {
            cnv.stotal = result.numberOfRanges
        }
        
        return result
    }

    private Regions extractSampleIndividualCnvs(String sample) {
        Regions sampleCNVs = new Regions()
        results.each { caller, calls ->
            for(Region cnv in calls.grep { it.sample == sample }) {
                sampleCNVs.addRegion(cnv)
            }
        }
        return sampleCNVs
    }
    
    /**
     * @param callers
     * @param cnv
     * @return true if the given CNV should be filtered out of the results
     */
    boolean isFiltered(List<String> callers, Region cnv) {
        
        // At least one of the CNV callers must pass quality filtering
        Map qualityFiltering = callers.collectEntries { caller ->
            [caller, checkCallQuality(cnv,caller)]    
        }
        
        if(!qualityFiltering.any { it.value == false }) {
            log.info "CNV $cnv filtered out by low quality score in $qualityFiltering"
            return true
        }
            
        // If a list of genes to filter to is set, include the CNV only if it
        // overlaps the gene list
        if(this.filterToGenes || this.excludeGenes) {
            List<String> cnvGenes = cnv.genes.tokenize(',')
            if(this.filterToGenes) {
                if(!cnvGenes.any { it in filterToGenes }) {
                    log.info "CNV $cnv filtered out because no genes overlap set gene list"
                    return true
                }
            }
            
            if(this.excludeGenes) {
                if(cnvGenes.any { it in excludeGenes}) {
                    log.info "CNV $cnv is filtered out because it overlaps an excluded gene"
                    return true
                }
            }
        }
        
        if(this.minimumCategory > 0) {
            if(!cnv.category)
                return true
            if(cnv.category<this.minimumCategory)
                return true
        }
        
        return false
    }
    
    /**
     * Return true if the caller failed quality check for the CNV
     * 
     * @param cnv
     * @param caller
     */
    def checkCallQuality(Region cnv, caller) {
        if(cnv[caller] == null) // Caller did not call this CNV
            return null

        if(qualityFilters[caller] == null)  // no threshold set for this caller
            return false

        if(cnv[caller].quality > qualityFilters[caller])
            return false

        return true // failed quality, below threshold
    }
    
    /**
     * Add contextual information to the CNV including:
     * <li>Type (del,dup)
     * <li>callers that found it
     * <li>number of callers that found it
     * <li>genes that the CNV overlaps
     * <li>variants that overlap the CNV / genes the CNV overlaps
     * @param sample
     * @param callers
     * @param cnv
     */
    void annotateCNV(String sample, List callers, Region cnv) {
        
        cnv.sample = sample
        
        List foundInCallers = []
        for(String caller in callers) {
            if(annotateCaller(cnv, caller)) {
                foundInCallers << caller
            }
        }
            
        List types = foundInCallers.collect { cnv[it]?.best?.type }.grep { it }.unique()
        if(types.size()>1)
            log.info "WARNING: CNV $cnv has conflicting calls: " + types
            
        cnv.type = types.join(",")
        cnv.count = callers.grep { it != "truth" }
                           .count { cnv[it].best != null } 
        
        annotateGenes(cnv)
            
        // Annotate the variants if we have a VCF for this sample
        if(variants[sample]) 
            sample = annotateVariants(sample, cnv)
            
        log.info "$cnv.count callers found $cnv of type $cnv.type in sample $sample covering genes $cnv.genes with ${cnv.variants?.size()} variants"
    }

    private void annotateGenes(Region cnv) {
        List<String> genes
        if(refGenes != null) {
            log.info "Annotating $cnv using RefGene database"
            
            if(!cnv.chr.startsWith('chr'))
                cnv = new Region('chr' + cnv.chr, cnv.range)
            
            genes = refGenes.getGenes(cnv)
            
            cnv.cdsOverlap = refGenes.getCDS(cnv)*.value?.sum()?:0
            
            log.info "CDS Overlap for $cnv is $cnv.cdsOverlap"
        }
        else {
            genes = targetRegions.getOverlaps(cnv)*.extra.unique()
        }

        
        cnv.genes=genes.join(",")

        annotateGeneCategories(cnv, genes)
    }
    
    /**
     * Annotates a CNV with how many calls support the CNV from different callers.
     * 
     * Support is divided into two levels: <code>supporting</code>, which means that there
     * is at least a minimum threshold of mutual overlap (default 50%), and <code>all</code>
     * which counts every overlapping call. In general, only supporting calls should be considered
     * for filtering and evidence that a call is real.
     * 
     * @param sample
     * @param cnv
     * @param caller
     * 
     * @return  true if the CNV was supported by at least 1 call from the given caller
     */
    @CompileStatic
    boolean annotateCaller(Region cnv, String caller) {
        
        final String sample = cnv['sample']
        
        // log.info "Find best CNV call for $caller"
        Iterable<Region> callerCalls = results[caller].grep { Region call -> call['sample'] == sample && call.overlaps(cnv) }
            
        Collection<Region> mutualOverlapCalls = callerCalls.grep { Region call -> call.mutualOverlap(cnv) > mergeOverlapThreshold }
            
        Region best = findBestCall(mutualOverlapCalls, caller)
        
        final Map props = [
            best : best,
            supporting : mutualOverlapCalls,
            all : callerCalls
        ]
            
        cnv[caller] = props        
            
        return !mutualOverlapCalls.isEmpty()
    }

    /**
     * Return the highest quality call from a list of CNVs 
     * 
     * @param mutualOverlapCalls
     * @param caller
     * @return
     */
    private Region findBestCall(Iterable<Region> mutualOverlapCalls, String caller) {
        final  Region best = mutualOverlapCalls.max { Region call -> call['quality']?.toFloat() }
        if(best != null) {
            log.info "Best CNV for $caller is " + best + " with quality " + best.quality
        }
        return best
    }

    private void annotateGeneCategories(Region cnv, List genes) {
        String sample = cnv.sample
        if(genelists) {
            String sampleGenelist = this.sampleToGenelist[sample]

            Map<String,Integer> categories = sampleGenelist ?
                    this.genelistCategories.get(sampleGenelist, this.geneCategories /* global */)
                    :
                    this.geneCategories

            cnv.category = genes.collect { categories[it] }.grep { it != null }.max()

            log.info "CNV $cnv assigned category $cnv.category based on $genes"
        } else {
            cnv.category = 0
        }
    }

    private String annotateVariants(String sample, Region cnv) {
        int sampleIndex = variants[sample][0].extra.header.samples.indexOf(sample)
        cnv.variants = variants[sample].getOverlaps(cnv)*.extra.grep { it.sampleDosage(sample) > 0 }
        cnv.hom = cnv.variants.count { it.dosages[sampleIndex] == 2 }
        cnv.het = cnv.variants.count { it.dosages[sampleIndex] == 1 }
        cnv.bal = Stats.mean(cnv.variants.grep { it.dosages[sampleIndex] == 1 }.collect {
            def refDepth = it.getAlleleDepths(0)[sampleIndex]
            if(refDepth == 0)
                return Float.NaN

            it.getAlleleDepths(1)[sampleIndex] / (float)refDepth
        })
        return sample
    }
}
