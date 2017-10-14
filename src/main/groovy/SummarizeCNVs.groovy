// vim: sw=4 expandtab cindent ts=4
import groovy.text.SimpleTemplateEngine
import groovy.util.logging.Log;

import org.codehaus.groovy.runtime.StackTraceUtils;
import org.omg.CORBA.SystemException

import graxxia.Matrix;

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
            truth 'Postiive control CNVs', args:1
            vcf 'VCF file containing variants for a sample in results', args:Cli.UNLIMITED
            target 'Target regions with id for each region to annotate', args:1, required:true
            chr 'Process CNVs only for given chromosome', args:1
            report 'Template to use for creating HTML report', args: 1
            dgv 'Path to file containing CNVs from database of genomic variants (DGV)', args: 1
            samples 'Samples to export', args:1
            quality 'Filtering by quality score, in form caller:quality...', args:Cli.UNLIMITED
            bam 'BAM file for a sample (if provided, used to customize IGV access)', args:Cli.UNLIMITED
            bampath 'Path to BAM files to enable access using IGV (can be URL)', args:1
            name 'Name for report (displayed in header)', args:1
            refgene 'Path to RefGene file downloaded from UCSC to annotate genes (optional)', args:1
            tsv 'Write consolidated CNV report in tab separated format to <file>', args:1
            imgpath 'Additional path to add for images', args:1
            genome 'Genome build to use when annotating CNVs', args:1
            idmask 'Mask to apply to sample ids for presentation in report', args:1
            genefilter 'Optional file of genes to filter CNVs to', args:1
            exgenes 'Optional file of genes to exclude from output', args:1
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
        try {
            
            RefGenes refGenes = RefGenes.download(opts.genome?:"hg19")
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
            
            if(opts.dgv) {
                summarizer.cnvAnnotator = new TargetedCNVAnnotator(target, opts.dgv)
            }
            
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
            
            Regions cnvs = summarizer.run(exportSamples)
            if(opts.tsv) {
                summarizer.writeTSV(cnvs, opts.tsv)
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
        
        if(this.dgvFile)
            this.cnvAnnotator = new TargetedCNVAnnotator(targetRegions, dgvFile)
        
        extractSamples()
        
        // Try to find a VCF file for each sample
        variants = samples.collectEntries { String sample ->
            [sample, vcfs.find { sample in it.samples }?.toRegions()]
        }

        if(exportSamples == null)
            exportSamples = this.samples

        Regions results = new Regions()
        for(s in exportSamples) {
            Regions sampleCnvs = extractSampleCnvs(s)
            for(cnv in sampleCnvs) {
                log.info "Merging CNV $cnv for sample '$s' type = $cnv.type"

                results.addRegion(cnv)
            }
        }
        
        // Second phase annotation (annotations that depend on the annotations 
        // created in 1st phase
        
        metaAnnotate(results)
        
        return results
    }
    
    void metaAnnotate(Regions results) {
        for(Region cnv in results) {
            cnv.samples = results.grep { it.overlaps(cnv) }*.sample.unique()
            cnv.sampleCount = cnv.samples.size()
            cnv.sampleFreq = (cnv.sampleCount / (double)samples.size())
            
            log.info "Samples for $cnv = $cnv.samples (count = $cnv.sampleCount)"
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
            
            w.println((["chr","start","end","sample","genes", "type","count","stotal","sampleCount","sampleFreq"] + 
                       (cnvAnnotator ? ["spanning","spanningFreq"] : []) +
                       cnvCallers + 
                       cnvCallers.collect { it+"_qual" }).join("\t"))
            
            cnvs.eachWithIndex { cnv, i ->
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
                ] + (cnvAnnotator ? [annotations[i].spanning.size(), annotations[i].spanningFreq] : []) +
                cnvCallers.collect { caller ->
                    cnv[caller] ? "TRUE" : "FALSE"
                }  + cnvCallers.collect { caller ->
                    cnv[caller] ? cnv[caller].quality : 0
                } 
                
                w.println line.join("\t")
            }
        }
    }
    
    void writeReport(Regions cnvs, 
                     String name, 
                     String fileName, 
                     String reportTemplate="cnv_report.html", 
                     boolean inlineJs=true, 
                     List bamFiles = [], 
                     def bamFilePath=false, 
                     String imgpath="") {
        
        log.info "Using report template: " + reportTemplate
        
        File outputFile = new File(fileName).absoluteFile
        log.info "Output path = " + outputFile.absolutePath
        
        SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
        String jsFileName = new File(reportTemplate).name.replaceAll('\\.html$','\\.js')
        InputStream templateStream 
        String jsCode = null
        
        File cnvReportFile = new File(reportTemplate)
        
        HTMLAssetSource assetSource
        File outputDir = outputFile.parentFile
        
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
                ) ;
                                
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
                    idMask: idMask
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
    
    Regions extractSampleCnvs(String sample) {
        
        // First accumulate CNVs from all the callers into one merged object
        Regions merged = new Regions()
        results.each { caller, calls ->
            for(Region cnv in calls.grep { it.sample == sample }) {
                merged.addRegion(cnv)
            }
        }
        
        Regions flattened = merged.reduce()
        
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
        
        List foundInCallers = []
        for(String caller in callers) {
            // log.info "Find best CNV call for $caller"
            println "results for $caller are " + results[caller]
            Region best = results[caller].grep { it.sample == sample && it.overlaps(cnv) }.max { it.quality?.toFloat() }
            if(best != null) {
                log.info "Best CNV for $caller is " + best + " with quality " + best?.quality
                foundInCallers << caller
            }
            cnv[caller] = best
        }
            
        def types = foundInCallers.collect { cnv[it]?.type }.grep { it }.unique()
        if(types.size()>1)
            log.info "WARNING: CNV $cnv has conflicting calls: " + types
            
            
        cnv.type = types.join(",")
        
        cnv.sample = sample
            
        cnv.count = callers.grep { it != "truth" }.count { cnv[it] != null } 
            
        if(refGenes != null) {
            log.info "Annotating $cnv using RefGene database"
            cnv.genes = refGenes.getGenes(cnv).join(",")
        }
        else {
            cnv.genes = targetRegions.getOverlaps(cnv)*.extra.unique().join(",")
        }
        
        // Annotate the variants if we have a VCF for this sample
        if(variants[sample]) {
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
        }
            
        log.info "$cnv.count callers found $cnv of type $cnv.type in sample $sample covering genes $cnv.genes with ${cnv.variants?.size()} variants"
    }
}
