import java.nio.file.FileSystems;
import java.nio.file.Files
import java.nio.file.Path
import java.util.regex.Matcher;

import graxxia.Matrix
import groovy.text.SimpleTemplateEngine
import groovy.transform.CompileStatic;
import groovy.util.logging.Log
import groovyx.gpars.GParsPool;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMRecord;

/**
 * The main entry point for the Ximmer simulation framework
 * 
 * @author simon
 */
@Log
class Ximmer {
    
    /**
     * Map of long name to short name for each caller. Long names 
     * are used in the configuration file because they are more 
     * descriptive. Short names appear in plots, result files and other places where
     * space is constrained.
     */
    Map<String,String> callerIdMap = [
        exome_depth : "ed", 
        cnmops : "cnmops",
        xhmm: "xhmm",
        conifer: "cfr"
    ]
    
    ConfigObject cfg = null
    
    File outputDirectory = null
    
    Map<String, SAM> bamFiles = [:]
    
    Pedigrees pedigrees = null
    
    Random random = null
    
    /**
     * If set, this seed will be used to initialise random number generation, for 
     * reproducible simulations
     */
    Integer seed = null
    
    Regions targetRegion = null
    
    Regions excludeRegions = null
    
    String runDirectoryPrefix = "run"
    
    File cacheDirectory = new File("cache")
    
    File dgvMergedFile
    
    File hg19RefGeneFile
    
    boolean enableSimulation = true
    
    int deletionsPerSample = -1
    
    List<String> callerIds = null
    
    Ximmer(ConfigObject cfg, String outputDirectory, boolean simulate) {
        this.outputDirectory = new File(outputDirectory)
        this.cfg = cfg
        this.random = this.seed != null ? new Random(this.seed)  : new Random()
        this.enableSimulation = simulate
        
        if(cfg.containsKey('deletionsPerSample'))
            this.deletionsPerSample = cfg.deletionsPerSample
        else
            this.deletionsPerSample = 1
            
        this.callerIds = ((Map<String,Object>)cfg.callers).keySet().collect { 
            if(callerIdMap[it] == null)
                throw new RuntimeException("Unknown CNV caller " + it + " referenced in configuration.")
            callerIdMap[it]
        }
    }
    
    List<File> runDirectories = []
    
   
    void run(analyse=true) {
        
        this.checkConfig()
        
        this.cacheReferenceData()
        
        this.resolveBamFiles()
        
        this.resolvePedigrees()
        
        this.targetRegion = new BED(cfg.target_regions).load().reduce()
        
        this.simulate()
        
        if(analyse) {
            List<String> analyses = this.runAnalysis()
        
            for(String analysis in analyses) {
                this.generateReport(analysis)
            }
        }
    }
    
    void checkConfig() {
        if(!cfg.containsKey('simulation_type')) 
            throw new RuntimeException("The key simulation_type is not found in the configuration file. Please set this to 'replace' or 'downsample'.")
            
        if(!(cfg.simulation_type in ["replace","downsample"])) 
            throw new RuntimeException("The key simulation_type is set to unknown value ${cfg.simulation_type}. Please set this to 'replace' or 'downsample'.")
    }
    
    void cacheReferenceData() {
        
        if(!cacheDirectory.exists())
            cacheDirectory.mkdirs()
            
        dgvMergedFile = new File(cacheDirectory, "dgvMerged.txt.gz")    
        if(!dgvMergedFile.exists()) {
            dgvMergedFile.withOutputStream { o -> 
                log.info("Downloading DGV database from UCSC ...")
                o << new URL("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz").openStream() 
            }
        }
        
        hg19RefGeneFile = new File(cacheDirectory, "refGene.txt.gz")    
        if(!hg19RefGeneFile.exists()) {
            hg19RefGeneFile.withOutputStream { o -> 
                log.info("Downloading DGV database from UCSC ...")
                o << new URL("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz").openStream() 
            }
        }
    }
    
    void simulate() {
        
        for(int i=0; i<cfg.runs; ++i) {
            String runDir = runDirectoryPrefix+i
            File dir = new File(outputDirectory, runDir)
            dir.mkdirs()
            runDirectories << dir
            if(!this.enableSimulation) {
                log.info "Simulation disabled: analysis will be performed directly from source files"
            }
            else {
                simulateRun(runDir)
            }
        }
    }
    
    List<String> runAnalysis() {
        List<String> analyses
        for(File runDir in runDirectories) {
            analyses = runAnalysisForRun(runDir)
        }
        return analyses
    }
    
    List<String> runAnalysisForRun(File runDir) {
        
        List<String> bamFiles 
        if(this.enableSimulation) {
            bamFiles = new File(runDir, "bams").listFiles().grep { File f -> f.name.endsWith(".bam") }.collect { "bams/"+it.name}
        } 
        else {
            bamFiles = this.bamFiles*.value*.samFile*.absolutePath
        }
        
        File bpipe = new File("eval/bpipe")
        
        String targetRegionsPath = new File(cfg.target_regions).absolutePath
        String toolsPath = new File("eval/pipeline/tools").absolutePath
        String ximmerSrc = new File("src/main/groovy").absolutePath
        int concurrency = cfg.containsKey("concurrency") ? cfg.concurrency : 2
        
        
        List<String> batches = ["analysis"]
        createAnalysis(runDir, batches[0])
        
        if(cfg.containsKey('analyses')) {
            for(String analysisName in cfg.analyses.keySet()) {
                // Create the corresponding analysis
                batches << createAnalysis(runDir, analysisName)
            }
        }
        
        if(batches.every { new File(runDir,it+"/report/cnv_report.html").exists()}) {
            log.info("Skipping bpipe run for $runDir because cnv_report.html already exists for all analyses (${batches.join(',')})")
            return batches
        }
        
        List<String> bpipeCommand = [
                "bash",
                bpipe.absolutePath,
                "run",
                "-n", "$concurrency",
                "-p","TOOLS=$toolsPath",
                "-p","DGV_CNVS=${dgvMergedFile.absolutePath}",
                "-p","XIMMER_SRC=$ximmerSrc",
                "-p", "callers=${callerIds.join(',')}",
                "-p", "refgene=${hg19RefGeneFile.absolutePath}",
                "-p", "simulation=${enableSimulation}",
                "-p", "batches=${batches.join(',')}",
                "-p", "target_bed=$targetRegionsPath", 
                "-p", "imgpath=${runDir.name}/analysis/report/",
                new File("eval/pipeline/exome_cnv_pipeline.groovy").absolutePath
                
            ] + bamFiles  + (enableSimulation ? ["true_cnvs.bed"] : [])
            
        log.info("Executing Bpipe command: " + bpipeCommand.join(" "))
        
        ProcessBuilder pb = new ProcessBuilder(
            bpipeCommand as String[]
        ).directory(runDir)
        
        Process p
        try {
            StringBuilder out = new StringBuilder()
            StringBuilder err = new StringBuilder()
            p = pb.start()
            
            p.waitForProcessOutput(out, err)
            int exitValue = p.waitFor()
            
            println out.toString()
            
            if(exitValue != 0) 
                throw new RuntimeException("Analysis failed for $runDir: " + err.toString())
            
        }
        finally {
          try { p.inputStream.close() } catch(Throwable t) { }      
          try { p.outputStream.close() } catch(Throwable t) { }      
          try { p.errorStream.close() } catch(Throwable t) { }      
        } 
        
        return batches
    }
    
    /**
     * Configure a bpipe analysis to run for the given named analysis
     * which must be a configured analysis either under 'callers' or
     * 'analyses' in the coniguration file. The analysis configured under
     * 'callers' is treated as a set of default values and these are customised
     * by entries under 'analyses'.
     * 
     * @param runDir
     * @param analysisName
     */
    String createAnalysis(File runDir, String analysisName) {
        
        // Clone the default configuration
        ConfigObject analysisCfg = cfg.callers.clone() 
                
        // If customised, merge the customised values for the analysis
        if(analysisName != "analysis") {
            analysisCfg.merge(cfg.analyses[analysisName])
            analysisName = "analysis-" + analysisName
        }
                
        
        File batchDir = new File(runDir, analysisName)
        batchDir.mkdirs()
        writeCallerParameterFile(batchDir, analysisCfg)
        return analysisName
    }
    
    /**
     * Write a file containing the caller parameters extracted from a config object
     * as the file caller.params.txt to the given directory.
     * 
     * @param outputDir
     * @param callersObj
     */
    void writeCallerParameterFile(File outputDir, ConfigObject callersObj) {
        String paramText = callersObj.collect { caller ->
            caller.value.collect { paramEntry ->
                "$paramEntry.key=$paramEntry.value"
            }

        }.flatten().join('\n')
        
        File paramFile = new File(outputDir, "caller.params.txt")
        paramFile.text = paramText
        log.info "Write file $paramFile with caller parameters"
    }
    
    void resolveBamFiles() {
        
        if(!cfg.containsKey("bam_files")) 
            throw new IllegalArgumentException("The configuration setting 'bam_files' is required")
            
        def bamFiles = cfg.bam_files
        if(bamFiles instanceof String) {
            bamFiles = MiscUtils.glob(bamFiles)
        }
        else 
        if(bamFiles instanceof List) {
            bamFiles = bamFiles.collect { MiscUtils.glob(it) }.flatten()
        }
        
        log.info("Resolved BAM files: " + bamFiles)
        
        this.bamFiles = bamFiles.collectEntries { bamPath ->
            SAM sam = new SAM(bamPath)
            
            if(sam.samples.size() > 1)
                throw new IllegalArgumentException("This program only supports single-sample BAM files")
                
            [sam.samples[0], new SAM(bamPath)]
        }
        
        if(cfg.containsKey("samples")) {
            Set<String> sampleSet =  (cfg.samples.males + cfg.samples.females) as Set
            
            this.bamFiles = this.bamFiles.grep { it.key in sampleSet }.collectEntries()
        }
        
    }
    
    void resolvePedigrees() {
        
        if(cfg.containsKey('ped_file')) {
            this.pedigrees = Pedigrees.parse(cfg.ped_file)
            return
        }
        
        if(!cfg.samples.containsKey('males') && !cfg.samples.containsKey('females'))
            throw new IllegalArgumentException("Please specify either a PED file in the configuration, or the samples.males and samples.femeales keys")
        
        this.pedigrees = new Pedigrees()
        if(cfg.samples.containsKey('males')) {
            Pedigrees males = Pedigrees.fromSingletons(cfg.samples.males)
            males.subjects*.value.each { ped -> ped.individuals.each { it.sex = Sex.MALE }}
            this.pedigrees.add(males)
        }
        
        if(cfg.samples.containsKey('females')) {
            Pedigrees females = Pedigrees.fromSingletons(cfg.samples.females)
            females.subjects*.value.each { ped -> ped.individuals.each { it.sex = Sex.FEMALE }}
            this.pedigrees.add(females)
        }
    }

    /**
     * Check that the sex of every sample can be determined
     */
    void checkPedigrees() {
        // TODO
    }
    
    void simulateByDownsample() {
    }
    
    List<Region> simulateRun(String runId) {
        
        // Make a directory for the run
        File runDir = new File(this.outputDirectory,  runId)
        File trueCnvsFile = new File(runDir,"true_cnvs.bed")
        
        if(trueCnvsFile.exists()) {
            log.info "Skipping generation of run $runId because $trueCnvsFile exists already!"
            return
        }
        
        runDir.mkdirs()
        if(!runDir.exists())
            throw new IOException("Unable to create run directory: $runDir")
         
        log.info "Created run directory " + runDir.absolutePath
        
        File bamDir = new File(runDir,"bams")
        bamDir.mkdirs()
        if(!bamDir.exists())
            throw new IOException("Unable to create run directory: $bamDir")
         
        int concurrency = cfg.containsKey("concurrency") ? cfg.concurrency : 2
        
        log.info "Creating output BAM files with concurrency $concurrency"
        
        // Check for existing bam files - by default we will not re-simulate for samples that
        // already have a bam file
        List<SAM> existingBams = Files.newDirectoryStream(bamDir.toPath(), '*.bam').grep { Path p ->
            // Ignore if the index for the file doesn't exist
            new File(p.toFile().absolutePath + '.bai').exists() 
        }.collect { Path p -> new SAM(p.toFile()) }
        
        List<String> existingSamples = Collections.synchronizedList(existingBams.collect { it.samples[0] })
        
        // Compute the target samples
        List<SAM> targetSamples = computeTargetSamples()
        List<Region> simulatedCNVs = []
        
        Ximmer me = this
        GParsPool.withPool(concurrency) {
            simulatedCNVs = targetSamples.collectParallel(me.&simulateSampleCNVs.curry(bamDir,existingSamples,existingBams))
        }
            
        trueCnvsFile.withWriter { w -> 
            w.println simulatedCNVs.collect { it as List }.flatten().collect { r -> 
                [r.chr, r.from, r.to+1, r.sample ].join("\t") 
            }.join("\n")
                
            if('known_cnvs' in cfg) {
                RangedData cnvs = new RangedData(cfg.known_cnvs).load(columnNames: ['chr','start','end','sample','type'])
                w.println cnvs.collect { r -> 
                    [r.chr, r.from, r.to+1, r.sample ].join("\t") 
                }.join("\n")
            }
        }
        return simulatedCNVs
    }
    
    Regions simulateSampleCNVs(File bamDir, List<String> existingSamples, List<SAM> existingBams, SAM targetSAM) {
        try {
            String sample = targetSAM.samples[0]
            Regions cnvs = checkExistingSimulatedSample(bamDir, existingSamples, existingBams, sample)
            if(cnvs != null) {
                log.info "Loaded pre-existing CNVs with samples " + cnvs*.sample
                return cnvs
            }
                          
            cnvs = simulateSampleCNV(bamDir, targetSAM)
            cnvs.each { it.sample = sample }
            return cnvs
        }
        catch(Exception e) {
            log.severe("Sample ${targetSAM.samples[0]} failed in simulation")
            log.throwing("Ximmer", "simulateRun", e)
            e.printStackTrace()
            throw e
        }
    }
    
    /**
     * Check if the given sample already has been simulated, and if it has,
     * return a Regions object containing the CNVs that were simulated. Otherwise
     * return null.
     * 
     * @param outputDir         the output directory where BAM files are placed
     * @param existingSamples   list of known samples extracted from BAMs - if not in there, will return null
     * @param existingBams      list of BAM files in output directory
     * @param sample            the sample to check for
     * 
     * @return      null if sample not simulated, otherwise Regions 
     *              object containing CNVs, with sample property set
     */
    Regions checkExistingSimulatedSample(File outputDir, List existingSamples, List existingBams, String sample) {
        
        int existingIndex = existingSamples.indexOf(sample)
        if(existingIndex<0) {
            return null
        }
        
        SAM existingSAM = existingBams[existingIndex]
        String sampleName = existingSAM.samFile.name
        
        File bedFile = new File(outputDir, sample + ".cnvs.bed")
        if(!bedFile.exists())
            return null
        
        log.info "Skipping simulation of CNV for sample $sample because a simulation BAM already exists for this sample in the output directory"
        
        Regions result = new Regions()
        
        Regions cnvs = new BED(bedFile).load()
        cnvs.each { 
            it.sample = sample 
            result.addRegion(it)
        }
        return result
    }
    
    Regions simulateSampleCNV(File outputDir, SAM targetSample) {
        
        String sampleId = targetSample.samples[0]
        Regions deletions = new Regions()
        Regions exclusions = this.excludeRegions ? this.excludeRegions.collect { it } as Regions : new Regions()
        
        for(int cnvIndex = 0; cnvIndex < this.deletionsPerSample; ++cnvIndex) {
            CNVSimulator simulator = new CNVSimulator(targetSample, null)
            if(this.seed != null) 
                simulator.random = new Random((this.seed<<16) + (cnvIndex << 8)  + (targetSample.hashCode() % 256))
  
            // Choose number of regions randomly in the range
            // the user has given
            int numRegions = cfg.regions.from + random.nextInt(cfg.regions.to - cfg.regions.from) 
            Region r = simulator.selectRegion(this.targetRegion, numRegions, exclusions)
            
            log.info "Seed region for ${sampleId} is $r" 
            
            exclusions.addRegion(r)
            deletions.addRegion(r)
        }
        
        File bedFile = new File(outputDir, sampleId + ".cnvs.bed")
        Region r = deletions[0]
        deletions.each { it.sample = sampleId }
        deletions.save(bedFile.absolutePath + ".tmp")
        
        File outputFile = new File(outputDir, sampleId + "_${r.chr}_${r.from}-${r.to}.bam" )
        
        CNVSimulator simulator = new CNVSimulator(targetSample, null)
        if(cfg.simulation_type == "downsample") {
                    
            simulator.simulationMode = "downsample"
        }
        else {
            throw new UnsupportedOperationException("Still working on this!")
        }
        
        simulator.createBam(outputFile.absolutePath, deletions)
        indexBAM(outputFile)
        
        deletions.save(bedFile.absolutePath)
        new File(bedFile.absolutePath + ".tmp").delete()
        return deletions
    }
    
    void indexBAM(File bamFile) {
            
        SAMFileReader reader = new SAMFileReader(bamFile)
        File outputFile = new File(bamFile.path + ".bai")
        BAMIndexer indexer = new BAMIndexer(outputFile, reader.getFileHeader());
            
        reader.enableFileSource(true);
        int totalRecords = 0;
            
        // create and write the content
        for (SAMRecord rec : reader) {
            indexer.processAlignment(rec);
        }
        indexer.finish();
    }
    
    /**
     * For an x-replacement simulation the target samples are the maximised set of 
     * unrelated females. For a downsample simulation, all the samples can be
     * targets.
     * 
     * @return
     */
    List<SAM> computeTargetSamples() {
        if(cfg.simulation_type == "downsample") {
            return this.bamFiles*.value
        }
        else {
            throw new UnsupportedOperationException("Still working on this!")
        }        
    }
    
    void generateReport(String analysisName) {
        
        log.info "Generating consolidated report for " + this.runDirectories.size() + " runs in analysis $analysisName"
        
        String summaryHTML = generateSummary(analysisName)
        
        generateROCPlots(analysisName)
        
        log.info("Generating HTML Report ...")
        File mainTemplate = new File("src/main/resources/index.html")
        
        String outputName = analysisName + ".html"
        
        new File(outputDirectory, outputName).withWriter { w ->
            SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
            templateEngine.createTemplate(mainTemplate.newReader()).make(
                analysisName : analysisName,
                runDirectories: runDirectories,
                outputDirectory : outputDirectory.name,
                summaryHTML : summaryHTML,
                callers: this.callerIds
            ).writeTo(w)
        }
        
        // We need to also copy the cnv.js file, because the individual CNV reports 
        // put it into the wrong location
        
        File cnvJs = new File(outputDirectory,"cnv.js")
        if(!cnvJs.exists()) {
            Files.copy(new File("src/main/resources/cnv_report.js").toPath(), 
                       cnvJs.toPath())
        }
    }
    
    
    /**
     * Return a list of simulated CNVs - each cnv has a 'sample' and 'run' 
     * property defining which run it was simulated in and which sample within 
     * that run.
     */
    List<Region> readCnvs() {
        List<Region> cnvs = []
        for(runDir in runDirectories) {
            File cnvFile = new File(runDir,"true_cnvs.bed")
            if(cnvFile.exists()) {
                new BED(cnvFile,withExtra:true).load().each { Region r ->
                    r.run = runDir.name
                    r.sample = r.extra
                    cnvs << r
                }
            }
        }               
        
        // Find the amount of each cnv intersected by the target regions
        cnvs.each { Region r ->
            // Number of targeted base pairs overlapped
            r.targeted = targetRegion.intersect(r)*.size().sum()
            
            // Number of targets overlapped
            r.targets = targetRegion.intersect(r).size()
        }
        
        return cnvs
    }
    
    File writeCombinedCNVInfo(List<Region> cnvs) {
        
       if(!outputDirectory.exists())
           outputDirectory.mkdirs()
       
       File cnvFile = new File(outputDirectory,"combined_cnvs.tsv")
       cnvFile.withWriter { w ->
            w.println([ "chr","start","end","sample","tgbp","targets" ].join("\t"))
            
            if(!enableSimulation) 
                return   
            
            cnvs.collect { 
                [
                    it.chr, 
                    it.from, 
                    it.to,
                    it.sample, 
                    it.targeted,
                    it.targets
                ]   
            }.each {
                w.println it.join("\t")
            }
        }
       return cnvFile
    }
    
    /**
     * Returns HTML for the summary / overview tab
     * 
     * @return
     */
    String generateSummary(String analysisName) {
        
        List<Region> cnvs = readCnvs()
        
        // Load the results
        List<RangedData> results = runDirectories.collect { runDir ->
            new RangedData(new File(runDir,"$analysisName/report/cnv_report.tsv").path).load([:], { r ->
                    callerIds.each { r[it] = r[it] == "TRUE" }
          })
        }
        
        plotCNVSizeHistograms(cnvs)
        
//        List<IntRange> sizeBins = [0..<200, 200..<500, 500..<1000, 1000..<2000,2000..<10000]
//        Map<IntRange,Integer> binnedCounts = 
        
        Map<String,Integer> callerCounts = this.callerIds.collectEntries { caller -> 
            [caller, results.sum { r -> r.count { cnv -> cnv.truth && cnv[caller] } }] 
        }
        
        File mainTemplate = new File("src/main/resources/summary.html")
        
        StringWriter result = new StringWriter()
        SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
        templateEngine.createTemplate(mainTemplate.newReader()).make(
            cnvs: cnvs,
            bamFiles: this.bamFiles,
            batch_name : outputDirectory.name,
            callers: this.callerIds,
            simulation_type: cfg.simulation_type,
            results: results,
            callerCounts: callerCounts
        ).writeTo(result)
        return result.toString()
    }
    
    void plotCNVSizeHistograms(List<Region> cnvs) {
        
        if(!this.enableSimulation)
            return
        
        File combinedCnvs = writeCombinedCNVInfo(cnvs)
        
        runPython([:], outputDirectory, new File("src/main/python/cnv_size_histogram.py"), [combinedCnvs.absolutePath, new File(outputDirectory,"cnv_size_histogram.png").absolutePath])
    }
    
    void generateROCPlots(String analysisName) {
        
       if(!this.enableSimulation)
            return
            
        log.info "Creating combined report for callers " + this.callerIds.join(",")
  
        runR(outputDirectory, 
            new File("src/main/R/ximmer_cnv_plots.R"), 
                ANALYSIS: analysisName,
                SRC: new File("src/main/R").absolutePath, 
                XIMMER_RUNS: runDirectories.size(),
                TARGET_REGION: new File(cfg.target_regions).absolutePath,
                XIMMER_CALLERS: this.callerIds.join(",")
            )
    }
    
    void runPython(Map env, File dir, File scriptFile, List<String> args=[]) {
        
        log.info("Executing python script ...")
        String rTempDir = MiscUtils.createTempDir().absolutePath
            
        String python = "python" 
        if(cfg.tools.containsKey("python"))
            python = cfg.tools.python
            
        File tempScript = new File(rTempDir, "XimmerPythonScript.sh")
        tempScript.setExecutable(true)
        
        tempScript.text = 
            """$python ${scriptFile.absolutePath} ${args.join(' ')}\n"""
            
        ProcessBuilder pb = new ProcessBuilder([
           "bash", tempScript.absolutePath
        ] as String[]).directory(dir)
        
        Map pbEnv = pb.environment()
        env.each { k,v ->
            pbEnv[k] = String.valueOf(v)
        }
                      
        Process process = pb.start()
        
        File pythonLogFile = new File("python.log")
        pythonLogFile.withOutputStream { o ->
            process.consumeProcessOutput(System.out, System.err)
        }
            
        int exitCode = process.waitFor()
        if(exitCode != 0) 
            throw new RuntimeException("Execution of python script $tempScript.absolutePath failed: \n"  + pythonLogFile.text + "\n")
    }
    
    
    void runR(Map env, File dir, File scriptFile) {
        
        log.info("Executing R script ...")
        String rTempDir = MiscUtils.createTempDir().absolutePath
            
        String rScriptExe = "Rscript" 
        if(cfg.tools.containsKey("Rscript"))
            rScriptExe = cfg.tools.Rscript
            
        File tempScript = new File(rTempDir, "XimmerRscript.R")
        tempScript.text = 
            """unset TMP; unset TEMP; TEMPDIR="$rTempDir" $rScriptExe - < ${scriptFile.absolutePath}\n"""
            
        tempScript.setExecutable(true)
        
        ProcessBuilder pb = new ProcessBuilder([
           "bash", tempScript.absolutePath
        ] as String[]).directory(dir)
        
        Map pbEnv = pb.environment()
        env.each { k,v ->
            pbEnv[k] = String.valueOf(v)
        }
                      
        Process process = pb.start()
        
        File rLogFile = new File("R.log")
        rLogFile.withOutputStream { o ->
            process.consumeProcessOutput(System.out, System.err)
        }
            
        int exitCode = process.waitFor()
        if(exitCode != 0) 
            throw new RuntimeException("Execution of R script $scriptFile failed: \n"  + rLogFile.text + "\n")
    }
    
    public static void main(String [] args) {
        
        println "*" * 50
        println """
        #     #  ###  #     #  #     #  #######  ######   
         #   #    #   ##   ##  ##   ##  #        #     #  
          # #     #   # # # #  # # # #  #        #     #  
           #      #   #  #  #  #  #  #  #####    ######   
          # #     #   #     #  #     #  #        #   #    
         #   #    #   #     #  #     #  #        #    #   
        #     #  ###  #     #  #     #  #######  #     #  
        """.stripIndent()
        println "*" * 50
                
        Cli cli = new Cli(usage:"ximmer -c <config> -o <output_directory>")
        cli.with {
            c "Configuration file", args:1, required:true
            o "Output directory", args:1, required:true
            seed "Random seed", args:1
            simonly "Run only simulation component"
            nosim "Analyse samples directly (no simulation)"
            v "Verbose output"
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        MiscUtils.configureSimpleLogging("ximmer.log")
        if(opts.v) {
            MiscUtils.configureVerboseLogging()
            println "Configured verbose logging"
            log.info "Configured verbose logging"
        }
        
        // Parse the configuration
        ConfigObject cfg = new ConfigSlurper().parse(new File(opts.c).text)
        
        Ximmer ximmer = new Ximmer(cfg, opts.o, !opts.nosim)
        if(opts.seed != false) {
            ximmer.seed = opts.seed.toInteger()
        }
        
        ximmer.run(!opts.simonly)
    }
}
