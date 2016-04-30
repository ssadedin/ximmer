import graxxia.Matrix
import groovy.text.SimpleTemplateEngine
import groovy.transform.CompileStatic;
import groovy.util.logging.Log
import net.sf.samtools.BAMIndexer;
import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMRecord;

/**
 * The main entry point for the Ximmer simulation framework
 * 
 * @author simon
 */
@Log
class Ximmer {
    
    ConfigObject cfg = null
    
    File outputDirectory = null
    
    Map<String, SAM> bamFiles = [:]
    
    Pedigrees pedigrees = null
    
    Random random = null
    
    Regions targetRegion = null
    
    Regions excludeRegions = null
    
    String runDirectoryPrefix = "run"
    
    File cacheDirectory = new File("cache")
    
    File dgvMergedFile
    
    File hg19RefGeneFile
    
    boolean enableSimulation = true
    
    Ximmer(ConfigObject cfg, String outputDirectory, boolean simulate) {
        this.outputDirectory = new File(outputDirectory)
        this.cfg = cfg
        this.random = new Random(42) // todo: programmable seed
        this.enableSimulation = simulate
    }
    
    List<File> runDirectories = []
    
   
    void run() {
        
        this.cacheReferenceData()
        
        this.resolveBamFiles()
        
        this.resolvePedigrees()
        
        this.targetRegion = new BED(cfg.target_regions).load().reduce()
        
        this.simulate()
        
        this.runAnalysis()
        
        this.generateReport()
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
    
    void runAnalysis() {
        for(File runDir in runDirectories) {
            runAnalysisForRun(runDir)
        }
    }
    
    void runAnalysisForRun(File runDir) {
        
        List<String> bamFiles 
        if(this.enableSimulation) {
            bamFiles = new File(runDir, "bams").listFiles().grep { File f -> f.name.endsWith(".bam") }.collect { "bams/"+it.name}
        } 
        else {
            bamFiles = this.bamFiles*.value*.samFile*.absolutePath
        }
        
        String callers = ((Map<String,Object>)cfg.callers).keySet().join(",")
        
        if(new File(runDir,"analysis/report/cnv_report.html").exists()) {
            log.info("Skipping bpipe run for $runDir because cnv_report.html already exists")
            return
        }
        
        File bpipe = new File("eval/bpipe")
        
        String targetRegionsPath = new File(cfg.target_regions).absolutePath
        String toolsPath = new File("eval/tools").absolutePath
        String ximmerSrc = new File("src/main/groovy").absolutePath
        
        List<String> bpipeCommand = [
                "bash",
                bpipe.absolutePath,
                "run",
                "-p","TOOLS=$toolsPath",
                "-p","DGV_CNVS=${dgvMergedFile.absolutePath}",
                "-p","XIMMER_SRC=$ximmerSrc",
                "-p", "callers=$callers",
                "-p", "refgene=${hg19RefGeneFile.absolutePath}",
                "-p", "simulation=${enableSimulation}",
                "-p", "batch_name=analysis",
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
         
        // Compute the target samples
        List<SAM> targetSamples = computeTargetSamples()
        List<Region> simulatedCNVs = []
        for(SAM targetSample in targetSamples) {
            // Choose number of regions randomly in the range
            // the user has given
            int numRegions = cfg.regions.from + random.nextInt(cfg.regions.to - cfg.regions.from)
            Region cnv = simulateSampleCNV(bamDir, targetSample, numRegions)
            cnv.sample = targetSample.samples[0]
            simulatedCNVs << cnv
        }
        
        trueCnvsFile.withWriter { w -> 
            w.println simulatedCNVs.collect { [it.chr, it.from, it.to+1, it.sample ].join("\t") }.join("\n")
        }
        return simulatedCNVs
    }
    
    Region simulateSampleCNV(File outputDir, SAM targetSample, int numRegions) {
        
        if(cfg.simulation_type == "downsample") {
            CNVSimulator simulator = new CNVSimulator(targetSample, null)
            
            String sampleId = targetSample.samples[0]
            Region r = simulator.selectRegion(this.targetRegion, numRegions, excludeRegions)
            log.info "Seed region for ${sampleId} is $r" 
                
            File outputFile = new File(outputDir, sampleId + "_${r.chr}_${r.from}-${r.to}.bam" )
            simulator.simulationMode = "downsample"
            simulator.createBam(outputFile.absolutePath, r)
            indexBAM(outputFile)
            r.sample = sampleId
            return r
        }
        else {
            throw new UnsupportedOperationException("Still working on this!")
        }
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
    
    void generateReport() {
        
        log.info "Generating consolidated report for " + this.runDirectories.size() + " runs"
        
        String summaryHTML = generateSummary()
        
        generateROCPlots()
        
        log.info("Generating HTML Report ...")
        File mainTemplate = new File("src/main/resources/index.html")
        
        List<String> callers = ((Map<String,Object>)cfg.callers).keySet().collect { 
           it == "mops" ? "cnmops" : it
        }
        
        new File(outputDirectory, "index.html").withWriter { w ->
            SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
            templateEngine.createTemplate(mainTemplate.newReader()).make(
                runDirectories: runDirectories,
                batch_name : outputDirectory.name,
                summaryHTML : summaryHTML,
                callers: callers
            ).writeTo(w)
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
    String generateSummary() {
        
        List<Region> cnvs = readCnvs()
        
        List<String> callers = ((Map<String,Object>)cfg.callers).keySet().collect { 
           it == "mops" ? "cnmops" : it
        }
        
        // Load the results
        List<RangedData> results = runDirectories.collect { runDir ->
            new RangedData(new File(runDir,"analysis/report/cnv_report.tsv").path).load([:], { r ->
                    callers.each { r[it] = r[it] == "TRUE" }
          })
        }
        
        plotCNVSizeHistograms(cnvs)
        
//        List<IntRange> sizeBins = [0..<200, 200..<500, 500..<1000, 1000..<2000,2000..<10000]
//        Map<IntRange,Integer> binnedCounts = 
        
        Map<String,Integer> callerCounts = callers.collectEntries { caller -> 
            [caller, results.sum { r -> r.count { cnv -> cnv.truth && cnv[caller] } }] 
        }
        
        File mainTemplate = new File("src/main/resources/summary.html")
        
        StringWriter result = new StringWriter()
        SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
        templateEngine.createTemplate(mainTemplate.newReader()).make(
            cnvs: cnvs,
            bamFiles: this.bamFiles,
            batch_name : outputDirectory.name,
            callers: callers,
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
    
    void generateROCPlots() {
        
       if(!this.enableSimulation)
            return
  
        runR(outputDirectory, 
            new File("src/main/R/ximmer_cnv_plots.R"), 
                SRC: new File("src/main/R").absolutePath, 
                XIMMER_RUNS: runDirectories.size(),
                TARGET_REGION: new File(cfg.target_regions).absolutePath
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
            throw new RuntimeException("Execution of R script $tempScript.absolutePath failed: \n"  + rLogFile.text + "\n")
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
            nosim "Analyse samples directly (no simulation)"
            v "Verbose output"
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        MiscUtils.configureSimpleLogging("ximmer.log")
        if(opts.v)
            MiscUtils.configureVerboseLogging()
        
        // Parse the configuration
        ConfigObject cfg = new ConfigSlurper().parse(new File(opts.c).text)
        
        new Ximmer(cfg, opts.o, !opts.nosim).run()
    }
}
