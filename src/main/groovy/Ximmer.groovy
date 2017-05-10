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
        conifer: "cfr",
        codex: "cdx"
    ]
    
    /**
     * Convert the expanded form of the configuration id into the simplified
     * form that is used in the analysis (eg: exomedepth => ed).
     * 
     * @param caller (possibly followed by underscore)
     * @return  caller with id replaced by compact form
     */
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
    
    List<String> males = null
    
    List<String> females = null
    
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
    
    String ximmerBase
    
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
        
        if(cfg.isSet('run_directory_prefix')) {
            this.runDirectoryPrefix = cfg.run_directory_prefix
            log.info "Setting run directory prefix to $runDirectoryPrefix"
        }
        else {
            log.info "No run directory configured: run directory = " + this.runDirectoryPrefix 
        }
        
        ximmerBase=System.properties['ximmer.base']
    }
    
    boolean enableTruePositives = false
   
    void run(analyse=true) {
        
        this.checkConfig()
        
        this.cacheReferenceData()
        
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
    
    List vcfFiles = []
    
    void checkConfig() {
        if(!cfg.containsKey('simulation_type')) 
            throw new RuntimeException("The key simulation_type is not found in the configuration file. Please set this to 'replace', 'downsample', or 'none'.")
            
        if(!(cfg.simulation_type in ["replace","downsample","none"])) 
            throw new RuntimeException("The key simulation_type is set to unknown value ${cfg.simulation_type}. Please set this to 'replace', 'downsample' or 'none'.")
            
//        if(!cfg.containsKey("runs") || !String.valueOf(cfg.runs).isInteger())
//            throw new RuntimeException("The key 'runs' must be set to an integer in the configuration file (current value is ${cfg.runs})")
            
        this.enableTruePositives = this.enableSimulation || ('known_cnvs' in cfg)
        
        if(cfg.containsKey('variants')) {
            if(cfg.variants instanceof List) {
                vcfFiles = cfg.variants
            }
            else
            if(cfg.variants instanceof String) {
                vcfFiles = cfg.variants.tokenize(',')*.trim()
            }
            else {
                throw new IllegalArgumentException("Configuration specificies variants using unsupported type " + cfg.variants.class.name)
            }
            
            List invalidVCFs = vcfFiles.grep { !new File(it).exists() }
            if(invalidVCFs)
                throw new IllegalArgumentException("The following VCF files were specified but do not exist: " + invalidVCFs.join(","))
        }
        
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
    
    Map<String,SimulationRun> runs
    
    void simulate() {
        
        this.runs = SimulationRun.configureRuns(this.outputDirectory, this.runDirectoryPrefix, cfg)
        
        this.bamFiles = this.runs.collect { it.value.bamFiles }.sum()
        
        for(Map.Entry<String,SimulationRun> runEntry in runs) {
            
            String runId = runEntry.key
            String trueCnvs = runEntry.value
            def runDir = runEntry.value.runDirectory
            if(!this.enableSimulation) {
                log.info "Simulation disabled: analysis will be performed directly from source files"
                File trueCnvsFile = new File(runDir,"true_cnvs.bed")
                trueCnvsFile.withWriter { w ->        
                    writeKnownCNVs(w, runId)
                }
            }
            else {
                simulateRun(runId)
            }
        }
    }
    
    List<AnalysisConfig> runAnalysis() {
        List<AnalysisConfig> analyses
        
        if(cfg.isSet('enable_parallel_analyses') && cfg.enable_parallel_analyses) {
            GParsPool.withPool(8) {      
                analyses = runs.values().collectParallel { run ->
                    runAnalysisForRun(run)
                }[0]
            }
            
        }
        else {
            for(SimulationRun run in runs.values()) {
                analyses = runAnalysisForRun(run)
            }
        }
        return analyses
    }
    
    Object analysisLock = new Object()
    
    /**
     * Execute Bpipe to perform an analysis for a simulation run.
     * 
     * @param run
     * @return
     */
    List<AnalysisConfig> runAnalysisForRun(SimulationRun run) {
        
        File runDir = run.runDirectory
        
        List<String> bamFiles 
        String targetRegionsPath 
        int concurrency
        List<AnalysisConfig> batches = []
        List drawCnvsParam = []
        
        synchronized(analysisLock) { // Avoid any potential multi-threading issues since all the below
                                     // are reading from non-threadsafe maps, config objects, etc.
            if(this.enableSimulation) {
                File bamDir = new File(run.runDirectory, "bams")
                log.info "Search for bam files in $bamDir"
                bamFiles = bamDir.listFiles().grep { File f -> f.name.endsWith(".bam") }.collect { "bams/"+it.name}
            } 
            else {
                bamFiles = run.bamFiles.collect { it.value.samFile.absolutePath }
            }
            
            assert !bamFiles.isEmpty() 
            
            targetRegionsPath = new File(cfg.target_regions).absolutePath
            concurrency = cfg.containsKey("concurrency") ? cfg.concurrency : 2 
           
            if(cfg.containsKey('analyses')) {
                for(String analysisName in cfg.analyses.keySet()) {
                    // Create the corresponding analysis
                    batches << createAnalysis(runDir, analysisName)
                }
            }
            else {
                batches << createAnalysis(runDir, 'analysis')
            }
            
            if(batches.every { new File(run.runDirectory,it.analysisName+"/report/cnv_report.html").exists()}) {
                log.info("Skipping bpipe run for $runDir because cnv_report.html already exists for all analyses (${batches*.analysisName.join(',')})")
                return batches
            }
            
            if(cfg.containsKey('draw_cnvs') && !cfg.draw_cnvs) {
                drawCnvsParam = ["-p","draw_cnvs=false"]
            }
        }
        
        
        String sampleIdMask = cfg.get('sample_id_mask','')
        
        File bpipe = new File("$ximmerBase/eval/bpipe")
        String toolsPath = new File("$ximmerBase/eval/pipeline/tools").absolutePath
        String ximmerSrc = new File("$ximmerBase/src/main/groovy").absolutePath
        
        List<String> bpipeCommand = [
                "bash",
                bpipe.absolutePath,
                "run",
                "-n", "$concurrency",
                "-p", "TOOLS=$toolsPath",
                "-p", "DGV_CNVS=${dgvMergedFile.absolutePath}",
                "-p", "XIMMER_SRC=$ximmerSrc",
                "-p", "callers=${callerIds.join(',')}",
                "-p", "refgene=${hg19RefGeneFile.absolutePath}",
                "-p", "simulation=${enableTruePositives}",
                "-p", "batches=${batches*.analysisName.join(',')}",
                "-p", "target_bed=$targetRegionsPath", 
                "-p", "sample_id_mask='$sampleIdMask'", 
                "-p", "imgpath=${runDir.name}/#batch#/report/", 
            ] + drawCnvsParam + [
                "$ximmerBase/eval/pipeline/exome_cnv_pipeline.groovy"
            ]  + bamFiles + vcfFiles + (enableTruePositives ? ["true_cnvs.bed"] : [])
            
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
        
        callerCfgs.addAll(writeCallerParameterFile(batchDir, cfg.callers, analysisCfg))
        
        List<String> callerLabels = callerCfgs.collect { caller ->
            if(analysisCfg[caller].containsKey('label')) {
                analysisCfg[caller].label
            }
            else {
                caller
            }
        }
        
        return new AnalysisConfig(analysisName:analysisName, callerCfgs: callerCfgs, callerLabels: callerLabels)
    }
    
    /**
     * Write a file containing the caller parameters extracted from a config object
     * as the file caller.params.txt to the given directory.
     * 
     * @param outputDir
     * @param defaultCfg        the default settings to apply to each caller
     * @param analysisConfig    the specific analysis to configure
     * 
     * @return a list of the analysis configs to run (caller ids with suffixes)
     */
    List<String> writeCallerParameterFile(File outputDir, ConfigObject defaultCfg, ConfigObject analysisConfig) {
        
       
        List<String> callerCfgs = []
        for(String key in analysisConfig.keySet()) {
            
            log.info "Analysis has caller configuration: $key"
            List callerParts = key.tokenize('_')
            String callerId = this.callerIdMap[callerParts[0]] ?:callerParts[0]
            if(callerParts.size() == 1)
                callerParts << "default"
            
            callerParts[0] = callerId
                
            ConfigObject callerParams = defaultCfg.clone()[callerId]
            
            log.info "Anaysis keys are: " + analysisConfig[key]*.key
            
            // Override defaults with analysis specific value
            analysisConfig[key].each { Map.Entry entry ->
                if(entry.key == "quality_filter") {
                    callerParams[key + "_" + entry.key] = entry.value
                }
                else {
                    callerParams[entry.key] = entry.value
                }
            }
            
            String paramText = callerParams.collect { paramEntry ->
                    "$paramEntry.key=$paramEntry.value"
            }.join('\n') + '\n'
            
            File paramFile = new File(outputDir, callerParts.join(".")+".params.txt")
            paramFile.text = paramText
            log.info "Write file $paramFile with caller parameters"
            
            callerCfgs << key
        }
        
        return callerCfgs
    }
    
    
    /**
     * Resolve pedigree / sample sex information.
     * <p>
     * Sample sex information can be specified in multiple different ways, so here
     * we check for all the possible configuration options and resolve everything
     * into a {@link pedigrees} attribute and {#males} and {#females} attributes that
     * are the authoritative source of that information.
     */
    void resolvePedigrees() {

        if(cfg.containsKey('ped_file')) {
            log.info "Reading pedigree information from file: " + cfg.ped_file
            this.pedigrees = Pedigrees.parse(cfg.ped_file)
            
            this.females = this.pedigrees.females
            this.males = this.pedigrees.males
        }
        else {
            if(!cfg.samples.containsKey('males') && !cfg.samples.containsKey('females')) {
                log.info "Sex of samples is not specified in configuration. Sex will be infered from data which can take additional processing time. Please consider adding sample sex to your configuration file"
                this.inferSexes()
            }
            createPedigreesFromSampleConfig()
        }
        
        // Initialize the final resolve set of males and females
        this.females = this.pedigrees.females
        this.males = this.pedigrees.males
    }
    
    /**
     * Convert pedigree information specified in the 'samples' configuration
     * block into an actual Pedigrees object.
     */
    void createPedigreesFromSampleConfig() {
        this.pedigrees = new Pedigrees()
        if(cfg.samples.containsKey('males')) {
            
            println "Males are " + cfg.samples.males
            
            Pedigrees males = Pedigrees.fromSingletons(cfg.samples.males)
            males.subjects*.value.each { ped -> ped.individuals.each { it.sex = Sex.MALE }}
            this.pedigrees = this.pedigrees.add(males)
        }
        
        if(cfg.samples.containsKey('females')) {
            Pedigrees females = Pedigrees.fromSingletons(cfg.samples.females)
            females.subjects*.value.each { ped -> ped.individuals.each { it.sex = Sex.FEMALE }}
            this.pedigrees = this.pedigrees.add(females)
        }
    }

    /**
     * Infer sexes bioinformatically from the data, and add them to the configuration.
     * <p>
     * If successful, the cfg.samples section will look as if the user specified the 
     * sexes themselves. Does not set the {@link #pedigrees} attribute.
     */
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
        File runDir = new File(this.outputDirectory,  this.runDirectoryPrefix + runId)
        
        log.info "Run directory = " + runDir
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
         
        log.info "Simulated BAMs will appear in $bamDir"
        
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
            knownCnvsFile = this.runs[runId].knownCnvs
        }
        else
        if('known_cnvs' in cfg) {
            knownCnvsFile = cfg.known_cnvs
        }
        
        if(!knownCnvsFile) {
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
        Regions simulationRegions = this.targetRegion
        if(cfg.simulation_type == "replace") {
            String maleId = this.males[random.nextInt(this.males.size())]
            sourceSample = this.bamFiles[maleId]
            if(sourceSample == null) 
                throw new IllegalStateException("The configured male sample $maleId does not have a corresponding BAM file configured under bam_files")
            simulationRegions = new Regions(this.targetRegion.grep { it.chr == 'chrX' || it.chr == 'X' })
        }
        
        if('groovy.lang.IntRange' != cfg.get('regions')?.class?.name)
            throw new IllegalArgumentException("The configured value for 'regions' is either missing or of incorrect type. Please set it to a range, for example: 1..5")
        
        for(int cnvIndex = 0; cnvIndex < this.deletionsPerSample; ++cnvIndex) {
            CNVSimulator simulator = new CNVSimulator(targetSample, sourceSample)
            if(this.seed != null) 
                simulator.random = new Random((this.seed<<16) + (cnvIndex << 8)  + (targetSample.hashCode() % 256))
  
            // Choose number of regions randomly in the range
            // the user has given
            int numRegions = cfg.regions.from + random.nextInt(cfg.regions.to - cfg.regions.from) 
            Region r = simulator.selectRegion(simulationRegions, numRegions, exclusions)
            
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
        
        assert !this.bamFiles.isEmpty()
        
        if(cfg.simulation_type == "downsample") {
            return this.bamFiles*.value
        }
        else {
            
            if(!("samples" in cfg) && !("ped_file" in cfg)) 
                throw new IllegalStateException("To use X replacement, please specify which samples are male and female in the 'samples' or 'ped_file' section of the configuration file")
            
            println "Females are " + this.pedigrees.females + " Males are " + this.pedigrees.males
            Set<String> females = this.pedigrees.females as Set
            
            List<SAM> results = this.bamFiles.grep { Map.Entry e ->
                println "Check: " + e.key + " is female in " + females
                e.key in females
            }*.value
        
            log.info "Females are " + females.join(',') + " with bam files " + this.bamFiles.keySet().join(',')
        
            if(results.isEmpty()) 
                throw new IllegalStateException("To use X replacement, there must be at least one female sample configured")
                
            return results
        }        
    }
    
    void generateReport(AnalysisConfig analysis) {
        
        // Find the analysable target region - this is actually produced by the bpipe run,
        // so it has to come from there
        
        File runDir = this.runs*.value[0].runDirectory
        File analysedTargets = new File(runDir,"target_regions.analysable.bed")
        
        log.info "Generating report based on analysed target regions: $analysedTargets"
        
        String analysisName = analysis.analysisName
        
        log.info "Generating consolidated report for " + this.runs.size() + " runs in analysis $analysisName"
        
        String summaryHTML = generateSummary(analysis)
        
        generateROCPlots(analysis, analysedTargets)
        
        log.info("Generating HTML Report ...")
        File mainTemplate = new File("$ximmerBase/src/main/resources/index.html")
        
        String outputName = analysisName + ".html"
        
        HTMLAssetSource source = new HTMLClassloaderAssetSource()
        HTMLAssets assets = new HTMLAssets(source, outputDirectory)
        
        assets << new HTMLAsset(source:'summary_report.js')
        assets << new HTMLAsset(source:'DOMBuilder.dom.min.js')
        
        String assetPayload = assets.render()
        
        new File(outputDirectory, outputName).withWriter { w ->
            SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
            templateEngine.createTemplate(mainTemplate.newReader()).make(
                analysisName : analysisName,
                runDirectories: runs*.value*.runDirectory,
                outputDirectory : outputDirectory.name,
                summaryHTML : summaryHTML,
                callers: this.callerIds,
                enableTruePositives: this.enableTruePositives,
                assets: assetPayload
            ).writeTo(w)
        }
        
        // We need to also copy the cnv.js file, because the individual CNV reports 
        // put it into the wrong location
        File cnvJs = new File(outputDirectory,"cnv.js")
        if(!cnvJs.exists()) {
            Files.copy(new File("$ximmerBase/src/main/resources/cnv_report.js").toPath(), 
                       cnvJs.toPath())
        }
        
        File cnvDiagramJs = new File(outputDirectory,"cnv_diagram.js")
        if(!cnvDiagramJs.exists()) {
            Files.copy(new File("$ximmerBase/src/main/resources/cnv_diagram.js").toPath(), 
                       cnvDiagramJs.toPath())
        } 
    }
    
    
    /**
     * Return a list of simulated CNVs - each cnv has a 'sample' and 'run' 
     * property defining which run it was simulated in and which sample within 
     * that run.
     */
    List<Region> readCnvs() {
        List<Region> cnvs = []
        for(File runDir in runs*.value*.runDirectory) {
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
        List<RangedData> results = runs*.value*.runDirectory.collect { runDir ->
            new RangedData(new File(runDir,"$analysisCfg.analysisName/report/cnv_report.tsv").path).load([:], { r ->
                    (callerIds + ["truth"]).each { r[it] = (r[it] == "TRUE") }
          })
        }
        
        plotCNVSizeHistograms(cnvs)
        
//        List<IntRange> sizeBins = [0..<200, 200..<500, 500..<1000, 1000..<2000,2000..<10000]
//        Map<IntRange,Integer> binnedCounts = 
        
        Map<String,Integer> callerCounts = this.callerIds.collectEntries { caller -> 
            [caller, results.sum { r -> r.count { cnv -> cnv.truth && cnv[caller] } }] 
        }
        
        File mainTemplate = new File("$ximmerBase/src/main/resources/summary.html")
        
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
        
        runPython([:], outputDirectory, new File("$ximmerBase/src/main/python/cnv_size_histogram.py"), 
                  [combinedCnvs.absolutePath, new File(outputDirectory,"cnv_size_histogram.png").absolutePath])
    }
    
    void generateROCPlots(AnalysisConfig analysisCfg, File analysedTargetRegions) {
        
       if(!enableTruePositives) {
           log.info "True positives not enabled: skipping ROC plots"
           return
       }
            
        String callerCfgs = analysisCfg.callerCfgs.collect { mapCallerId(it) }.unique().join(',')
        
        log.info "Creating combined ROC plots for callers " + callerCfgs + " with configurations " + analysisCfg.callerCfgs.join(",")
        
        runR(outputDirectory, 
            new File("$ximmerBase/src/main/R/ximmer_cnv_plots.R"), 
                ANALYSIS: analysisCfg.analysisName,
                SRC: new File("$ximmerBase/src/main/R").absolutePath, 
                XIMMER_RUNS: runs*.value*.runDirectory*.name.join(","),
                TARGET_REGION: analysedTargetRegions.absolutePath,
                XIMMER_CALLERS: callerCfgs,
                XIMMER_CALLER_LABELS: analysisCfg.callerLabels.join(','),
                SIMULATION_TYPE: cfg.simulation_type
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
        
        log.info "R environment: $env"
        
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
        
        if(opts.nosim && !cfg.containsKey('simulation_type'))
            cfg.simulation_type = 'none'
        
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
    
    List<String> callerLabels
    
    List<String> getCallerIds() {
        return callerCfgs.collect { String cfg -> cfg.tokenize('_')[0] }
    }
}


