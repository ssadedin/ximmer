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
        exomedepth : "ed", 
        cnmops : "cnmops",
        xhmm: "xhmm",
        conifer: "cfr"
    ]
    
    String mapCallerId(String caller) {
        for(String callerId in callerIdMap.keySet()) {
            if(caller.startsWith(callerId)) {
                return caller.replaceAll('^'+callerId, callerIdMap[callerId])
            }
        }
        return caller
    }
    
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
        this.enableSimulation = simulate && cfg.simulation_type != 'none'

        if(!enableSimulation)
            log.info "Simulation disabled!"
        
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
    
    boolean enableTruePositives = false
   
    void run(analyse=true) {
        
        this.checkConfig()
        
        this.cacheReferenceData()
        
        this.resolveBamFiles()
        
        this.targetRegion = new BED(cfg.target_regions).load().reduce()

        this.resolvePedigrees()
       
        this.simulate()
        
        if(analyse) {
            List<AnalysisConfig> analyses = this.runAnalysis()
        
            for(AnalysisConfig analysis in analyses) {
                this.generateReport(analysis)
            }
        }
    }
    
    void checkConfig() {
        if(!cfg.containsKey('simulation_type')) 
            throw new RuntimeException("The key simulation_type is not found in the configuration file. Please set this to 'replace', 'downsample', or 'none'.")
            
        if(!(cfg.simulation_type in ["replace","downsample","none"])) 
            throw new RuntimeException("The key simulation_type is set to unknown value ${cfg.simulation_type}. Please set this to 'replace', 'downsample' or 'none'.")
            
//        if(!cfg.containsKey("runs") || !String.valueOf(cfg.runs).isInteger())
//            throw new RuntimeException("The key 'runs' must be set to an integer in the configuration file (current value is ${cfg.runs})")
            
        this.enableTruePositives = this.enableSimulation || ('known_cnvs' in cfg)
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
    
    Map<String,String> runs
    
    
    void parseRuns() {
        // Map of run id to true_cnvs file
        def cfgRuns = cfg.runs
        if(cfgRuns instanceof String || cfgRuns instanceof Integer) {
            runs = (0..String.valueOf(cfgRuns).toInteger()).collectEntries { runId ->
                if('known_cnvs' in cfg)
                    [runId, cfg.known_cnvs]
                else
                    [runId, null]
            }
        }
        else
        if(cfgRuns instanceof ConfigObject) {
            log.info "Found named runs configured"
            runs = cfg.runs.collectEntries { [it.key, it.value.known_cnvs] }
            
            List missingRuns = runs.grep { !new File(it.value).exists() } 
            if(missingRuns)
                throw new Exception("One or more of the configured known_cnvs files does not exist: " + missingRuns*.value)
            
        }
        else 
            throw new Exception("The 'runs' configuration parameter is not in a recognised format. Expect either an integer or a child configuration, got $cfgRuns.")
    }
    
    void simulate() {
        
        this.parseRuns()
        
        for(Map.Entry runEntry in runs) {
            
            String runId = runEntry.key
            String trueCnvs = runEntry.value
            
            String runDir = runDirectoryPrefix+runId
            File dir = new File(outputDirectory, runDir)
            dir.mkdirs()
            runDirectories << dir
            if(!this.enableSimulation) {
                log.info "Simulation disabled: analysis will be performed directly from source files"
                File trueCnvsFile = new File(dir,"true_cnvs.bed")
                trueCnvsFile.withWriter { w ->        
                    writeKnownCNVs(w, runId)
                }
            }
            else {
                simulateRun(runDir)
            }
        }
    }
    
    List<AnalysisConfig> runAnalysis() {
        List<AnalysisConfig> analyses
        for(File runDir in runDirectories) {
            analyses = runAnalysisForRun(runDir)
        }
        return analyses
    }
    
    List<AnalysisConfig> runAnalysisForRun(File runDir) {
        
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
        
        
        List<AnalysisConfig> batches = []
//        batches << createAnalysis(runDir, "analysis") // default analysis
        
        if(cfg.containsKey('analyses')) {
            for(String analysisName in cfg.analyses.keySet()) {
                // Create the corresponding analysis
                batches << createAnalysis(runDir, analysisName)
            }
        }
        
        if(batches.every { new File(runDir,it.analysisName+"/report/cnv_report.html").exists()}) {
            log.info("Skipping bpipe run for $runDir because cnv_report.html already exists for all analyses (${batches*.analysisName.join(',')})")
            return batches
        }
        
        List drawCnvsParam = []
        if(cfg.containsKey('draw_cnvs') && !cfg.draw_cnvs) {
            drawCnvsParam = ["-p","draw_cnvs=false"]
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
                "-p", "simulation=${enableTruePositives}",
                "-p", "batches=${batches*.analysisName.join(',')}",
                "-p", "target_bed=$targetRegionsPath", 
                "-p", "imgpath=${runDir.name}/#batch#/report/", 
            ] + drawCnvsParam + [
                new File("eval/pipeline/exome_cnv_pipeline.groovy").absolutePath
            ]  + bamFiles  + (enableTruePositives ? ["true_cnvs.bed"] : [])
            
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
    AnalysisConfig createAnalysis(File runDir, String analysisName) {
        
        // Clone the default configuration
        ConfigObject analysisCfg = cfg.callers.clone() 
                
        if(analysisName != "analysis") {
            analysisCfg = cfg.analyses[analysisName]
            analysisName = "analysis-" + analysisName
        }
                
        
        File batchDir = new File(runDir, analysisName)
        batchDir.mkdirs()
        
        // The total list of all CNV caller configurations. 
        // A single caller might have multiple configurations within an analysis,
        // eg: xhmm_1, xhmm_2 etc.
        List<String> callerCfgs = []
        
        for(String caller in callerIds) {
            callerCfgs.addAll(writeCallerParameterFile(caller, batchDir, cfg.callers, analysisCfg))
        }
        
        return new AnalysisConfig(analysisName:analysisName, callerCfgs: callerCfgs)
    }
    
    /**
     * Write a file containing the caller parameters extracted from a config object
     * as the file caller.params.txt to the given directory.
     * 
     * @param outputDir
     * @param callersObj
     * @return a list of the analysis configs to run (caller ids with suffixes)
     */
    List<String> writeCallerParameterFile(String caller, File outputDir, ConfigObject defaultCfg, ConfigObject analysisConfig) {
        
       
        List<String> callerCfgs = []
        for(String key in analysisConfig.keySet()) {
            
            log.info "Analysis has caller configuration: $key"
            List callerParts = key.tokenize('_')
            String callerId = this.callerIdMap[callerParts[0]] ?:callerParts[0]
            if(callerParts.size() == 1)
                callerParts << "default"
            
            callerParts[0] = callerId
                
            ConfigObject callerParams = defaultCfg.clone()[callerId]
            
            // Override defaults with analysis specific value
            analysisConfig[key].each {
                callerParams[it.key] = it.value
            }
            
            String paramText = callerParams.collect { paramEntry ->
                    "$paramEntry.key=$paramEntry.value"
            }.join('\n')
            
            File paramFile = new File(outputDir, callerParts.join(".")+".params.txt")
            paramFile.text = paramText
            log.info "Write file $paramFile with caller parameters"
            
            callerCfgs << key
        }
        
        return callerCfgs
    }
    
    void addBamFiles(List<String> bamFilePaths) {
       
       Set knownSamples = (this.bamFiles*.key) as Set
        
       this.bamFiles += bamFilePaths.collectEntries { bamPath ->
            SAM sam = new SAM(bamPath)
            if(sam.samples.toUnique().size() > 1)
                throw new IllegalArgumentException("BAM file $bamPath contains mulitiple samples (${sam.samples.join(",")}). This program only supports single-sample BAM files")
            
            String sampleId = sam.samples[0]
            if(sampleId in knownSamples) 
                throw new IllegalArgumentException("Sample $sampleId appears in more than one BAM file. This program only supports a sample appearing in a single BAM file")
                
            knownSamples << sampleId
            
            [sampleId, new SAM(bamPath)]
        } 
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
        
        this.addBamFiles(bamFiles)
        
        if(cfg.containsKey("samples")) {
            Set<String> sampleSet =  (cfg.samples.males + cfg.samples.females) as Set
            
            this.bamFiles = this.bamFiles.grep { it.key in sampleSet }.collectEntries()
        }

        if(this.bamFiles.isEmpty())
            throw new RuntimeException("After comparing the samples specified in the 'samples' block with the available BAM files, no samples remain to analyse. Please check that sample ids are consistent with those in the BAM files supplied")
        
    }
    
    void resolvePedigrees() {

        if(cfg.containsKey('ped_file')) {
            this.pedigrees = Pedigrees.parse(cfg.ped_file)
            return
        }

        if(!cfg.samples.containsKey('males') && !cfg.samples.containsKey('females')) {
            log.info "Sex of samples is not specified in configuration. Sex will be infered from data which can take additional processing time. Please consider adding sample sex to your configuration file"
            this.inferSexes()
        }
        //    throw new IllegalArgumentException("Please specify either a PED file in the configuration, or the samples.males and samples.females keys")
        
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

    void inferSexes() {
        cfg.samples.females = []
        cfg.samples.males = []

        List<SexKaryotyper> unknownResults = []

        for(SAM sam in this.bamFiles*.value) {
            log.info "Checking sex of ${sam.samples[0]} ..."
            SexKaryotyper typer = new SexKaryotyper(sam, this.targetRegion)
            typer.run()
            if(typer.sex == Sex.FEMALE)
                cfg.samples.females << sam.samples[0]
            else
            if(typer.sex == Sex.MALE)
                cfg.samples.males << sam.samples[0]
            else
                unknownResults << typer
        }

        if(unknownResults) {
            throw new RuntimeException("""
                Unable to infer sex for one or more samples: 

                ${unknownResults.collect{it.bam.samples[0]+'('+it.sex.name()+')'}.join(', ')}

                Please add sexes manually in the samples block of the configuration file.
            """.stripIndent())
        }
    }

    /**
     * Check that the sex of every sample can be determined
     */
    void checkPedigrees() {
        // TODO
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
            
        writeCNVFile(trueCnvsFile, simulatedCNVs, runId)
        
        return simulatedCNVs
    }
    
    /**
     * Write out a file in BED-like format including all known true positive CNVs for this run.
     * These include both simulated CNVs and also any true positives known to be in the data.
     * 
     * @param trueCnvsFile
     * @param simulatedCNVs
     */
    void writeCNVFile(File trueCnvsFile, List<Region> simulatedCNVs, String runId) {
        trueCnvsFile.withWriter { w -> 
            w.println simulatedCNVs.collect { it as List }.flatten().collect { r -> 
                [r.chr, r.from, r.to+1, r.sample ].join("\t") 
            }.join("\n")
                
            writeKnownCNVs(w, runId)
        }
    }
    
    void writeKnownCNVs(Writer w, String runId) {
        
        String knownCnvsFile
        
        if(this.runs[runId]) {
            knownCnvsFile = this.runs[runId]
        }
        else
        if('known_cnvs' in cfg) {
            knownCnvsFile = cfg.known_cnvs
        }
        else {
            log.info "No true positive CNVs configured for run $runId"
            return
        }
        
        RangedData cnvs = new RangedData(knownCnvsFile).load(columnNames: ['chr','start','end','sample','type'])
        List knownCnvs = cnvs.collect { r -> 
            [r.chr, r.from, r.to+1, r.sample ].join("\t") 
        }
        w.println knownCnvs.join("\n")
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
        
        // If simulation mode is replacement, we need to select a male to simulate from
        SAM sourceSample = null
        if(cfg.simulation_type == "replace") {
            String maleId = cfg.samples.males[random.nextInt(cfg.samples.males.size())]
            sourceSample = this.bamFiles[maleId]
            if(sourceSample == null) 
                throw new IllegalStateException("The configured male sample $maleId does not have a corresponding BAM file configured under bam_files")
        }
        
        for(int cnvIndex = 0; cnvIndex < this.deletionsPerSample; ++cnvIndex) {
            CNVSimulator simulator = new CNVSimulator(targetSample, sourceSample)
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
        
        CNVSimulator simulator = new CNVSimulator(targetSample, sourceSample) 
        if(cfg.simulation_type == "downsample") {
           simulator.simulationMode = "downsample"
        }
        else {
           simulator.simulationMode = "replace"
           
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
            if(!("samples" in cfg)) 
                throw new IllegalStateException("To use X replacement, please specify which samples are male and female in the 'samples' section of the configuration file")
            
            Set<String> females = cfg.samples.females
            
            List<SAM> results = this.bamFiles.grep { Map.Entry e ->
                e.key in cfg.samples.females
            }*.value
        
            if(results.isEmpty()) 
                throw new IllegalStateException("To use X replacement, there must be at least one female sample configured")
                
            return results
        }        
    }
    
    void generateReport(AnalysisConfig analysis) {
        
        String analysisName = analysis.analysisName
        
        log.info "Generating consolidated report for " + this.runDirectories.size() + " runs in analysis $analysisName"
        
        String summaryHTML = generateSummary(analysis)
        
        generateROCPlots(analysis)
        
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
                callers: this.callerIds,
                enableTruePositives: this.enableTruePositives
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
            
            if(!enableTruePositives) 
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
    String generateSummary(AnalysisConfig analysisCfg) {
        
        List<Region> cnvs = readCnvs()
        
        // Load the results
        List<RangedData> results = runDirectories.collect { runDir ->
            new RangedData(new File(runDir,"$analysisCfg.analysisName/report/cnv_report.tsv").path).load([:], { r ->
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
        
        if(!enableTruePositives)
            return
        
        File combinedCnvs = writeCombinedCNVInfo(cnvs)
        
        runPython([:], outputDirectory, new File("src/main/python/cnv_size_histogram.py"), [combinedCnvs.absolutePath, new File(outputDirectory,"cnv_size_histogram.png").absolutePath])
    }
    
    void generateROCPlots(AnalysisConfig analysisCfg) {
        
       if(!enableTruePositives)
            return
            
        String callerCfgs = analysisCfg.callerCfgs.collect { mapCallerId(it) }.unique().join(',')
        
        log.info "Creating combined report for callers " + callerCfgs + " with configurations " + analysisCfg.callerCfgs.join(",")
  
        
        runR(outputDirectory, 
            new File("src/main/R/ximmer_cnv_plots.R"), 
                ANALYSIS: analysisCfg.analysisName,
                SRC: new File("src/main/R").absolutePath, 
                XIMMER_RUNS: runDirectories.join(","),
                TARGET_REGION: new File(cfg.target_regions).absolutePath,
                XIMMER_CALLERS: callerCfgs
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
        
        def additionalBams = opts.arguments()
        
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
        
        if(!additionalBams.isEmpty())
            ximmer.addBamFiles(additionalBams)
        
        ximmer.run(!opts.simonly)
    }
}

class AnalysisConfig {
    
    String analysisName
    
    List<String> callerCfgs
    
    List<String> getCallerIds() {
        return callerCfgs.collect { String cfg -> cfg.tokenize('_')[0] }
    }
}


