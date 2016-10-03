import groovy.util.logging.Log;

@Log
class SimulationRun {
    
    String id
    
    Map<String, SAM> bamFiles = [:]
    
    String knownCnvs
    
    File runDirectory
    
    void addBamFiles(List<String> bamFilePaths) {
       
       Set knownSamples = (this.bamFiles*.key) as Set
       
       log.info "Processing BAM Files: \n" + bamFilePaths.join('\n')
        
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
    
    void resolveBamFiles(def bamFileSpec, ConfigObject cfg) {
        
        def bamFiles = bamFileSpec
        
        if(bamFiles instanceof String) {
            bamFiles = MiscUtils.glob(bamFiles)
        }
        else 
        if(bamFiles instanceof List) {
            bamFiles = bamFiles.collect { MiscUtils.glob(it) }.flatten()
        }
        
        log.info("Resolved BAM files for run $id: " + bamFiles)
        
        this.addBamFiles(bamFiles)
        
        if(cfg.containsKey("samples")) {
            Set<String> sampleSet =  (cfg.samples.males + cfg.samples.females) as Set
            
            this.bamFiles = this.bamFiles.grep { it.key in sampleSet }.collectEntries()
        }

        if(this.bamFiles.isEmpty())
            throw new RuntimeException("After comparing the samples specified in the 'samples' block with the available BAM files, no samples remain to analyse. Please check that sample ids are consistent with those in the BAM files supplied")
        
    }
    
    
    static Map<String,SimulationRun> configureRuns(File outputDirectory, String runDirectoryPrefix, ConfigObject cfg) {
        
        // Map of run id to true_cnvs file
        def cfgRuns = 'runs' in cfg ? cfg.runs : false
        if(cfgRuns instanceof String || cfgRuns instanceof Integer) {
            configureIntegerRuns(outputDirectory, runDirectoryPrefix, cfgRuns, cfg)
        }
        else
        if(cfgRuns instanceof ConfigObject) {
            configureComplexRuns(outputDirectory, runDirectoryPrefix, cfg)
        }
        else 
            throw new Exception("The 'runs' configuration parameter is not in a recognised format. Expect either an integer or a child configuration, got $cfgRuns.")
    }

    private static Map<String,SimulationRun> configureComplexRuns(File outputDirectory, String runDirectoryPrefix, ConfigObject cfg) {
        
        Map<String,SimulationRun> runs = [:]
        
        String defaultKnownCnvs = cfg.isSet('known_cnvs') ? cfg.known_cnvs : null
        
        log.info "Found named runs configured"
        runs = cfg.runs.collectEntries { runEntry ->
            
            String runDir = runDirectoryPrefix + runEntry.key
            File dir = new File(outputDirectory, runDir)
            dir.mkdirs()
        
            SimulationRun run = 
                new SimulationRun(id:runEntry.key,
                                  knownCnvs:runEntry.value.isSet('known_cnvs') ? runEntry.value.known_cnvs : defaultKnownCnvs,
                                  runDirectory: dir)
                
            // Bam file could be specified globally (same for all runs, from cfg.bam_files)
            // or it could be specified specifically for this run
            def bamFileSpec  = runEntry.value.isSet('bam_files') ? runEntry.value.bam_files : cfg.bam_files
            
            if(bamFileSpec instanceof ConfigObject)
                throw new Exception("The 'bam_files' configuration element was not set. Please set this to a glob style path matching the bam files you wish to include")
            
            run.resolveBamFiles(bamFileSpec, cfg)
            
            [runEntry.key, run] 
        }
        
        log.info "Configured runs: " + runs.keySet().join(',')

        List missingRuns = runs.grep { it.value.knownCnvs && !new File(it.value.knownCnvs).exists() }
        if(missingRuns)
            throw new Exception("One or more of the configured known_cnvs files does not exist: " + missingRuns*.value*.knownCnvs)
            
        return runs
    }

    private static Map<String,SimulationRun> configureIntegerRuns(File outputDirectory, String runDirectoryPrefix, /* Integer or String */ Object cfgRuns, ConfigObject cfg) {
        
        String defaultKnownCnvs = cfg.isSet('known_cnvs') ? cfg.known_cnvs : null
        
        Map<String,SimulationRun> runs = [:]
        if(!String.valueOf(cfgRuns).isInteger())
            throw new IllegalArgumentException("The runs parameter (" + String.valueOf(cfgRuns) + ") must be an integer if not specified as a nested configuration.")
            
        if(cfg.bam_files instanceof ConfigObject)
            throw new Exception("The 'bam_files' configuration element was not set. Please set this to a glob style path matching the bam files you wish to include")
  
        
        runs = (0..String.valueOf(cfgRuns).toInteger()).collectEntries { runId ->
            
            String runDir = runDirectoryPrefix+runId
            File dir = new File(outputDirectory, runDir)
            dir.mkdirs()
             
            SimulationRun run = null
            if('known_cnvs' in cfg) {
                run = new SimulationRun(id:runId, knownCnvs:defaultKnownCnvs)
            }
            else {
                run = new SimulationRun(id:runId, knownCnvs:null)
            }
            
            // When configured as a simple integer, the bam files are assumed to be the same for all runs
            run.resolveBamFiles(cfg.bam_files, cfg)
            return [runId, run]
        }
        return runs
    }
  
}
