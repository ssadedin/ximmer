println "=" * 100
println "CNV Simulator " + (new Date().toString())
println "=" * 100
println ""

Cli cli = new Cli(usage:"SimulateCNVs [options]\n")
cli.with {
    n "Number of regions the CNVs should cover", args:1, required: true
    f "Sample names of females to simulate for (required if not using -ped)", longOpt: "females", args: 1
    m "Sample names of males from which to select from when simulating (required if not using -ped)", longOpt:"male", args: 1
    ped "PED file describing sexes of samples supplied by -bam (required unless using -f, -m)", args:1
    r "BED file of regions in which CNVs should be simulated", longOpt: "regions", args:1, required: true
    o "output file to write regions of true CNVs to (BED format)", longOpt: "output", args:1, required:true
	bam "BAM files for samples (optional if bam files are specified in sample information file)", args:1
    so "output file to write sample information to", args:1
    i "sample information file (tab separated, MGHA format, required for -so)", args:1
    cov "mean target coverage for simulated samples (default: maximise) - requires -abed", args:1
    abed "bed file of covered regions (required for -cov)", args:1
    dgv "File of population variants to exclude (UCSC Database of Genomic Variants format) (optional)", args:1
    dgvFreq "Population frequency threshold at which to exclude regions from simulation (0.01)", args:1
    d "directory to store output files in", longOpt: "dir", args:1
    seed "Random seed to use (for reproducible simulations)", args:1
	threads "Number of threads (concurrency) to use in simulating", args:1
    mode "Method for simulating deletions: downsampling reads=downsample, replacement from male=replace, both=use both methods (2 output files)", args:1
}

opts = cli.parse(args)
if(!opts)
    System.exit(1)

File outputDir = opts.d?new File(opts.d):null    
if(opts.d && !outputDir.exists())
    if(!outputDir.mkdirs() && !outputDir.exists()) // second check necessary due to race condition
        throw new IOException("Cannot create directory ${opts.d}")

def males = opts.ms
def females = opts.fs

if(!males || !females) {
    
    if(!opts.ped) {
        System.err.println "Please specify sexes of samples using either -f and -m options, or alternatively, using -ped"
        cli.usage()
        System.exit(1)
    }
    
    // Get sexes from PED file
    Pedigrees peds = Pedigrees.parse(opts.ped)
    males = peds.subjects.grep { it.value.sex == Sex.MALE }.collect { it.id }
    females = peds.subjects.grep { it.value.sex == Sex.FEMALE }.collect { it.id }
}

Map samples = [:]
if(opts.i)
  samples = SampleInfo.parse_sample_info(opts.i)

  
Regions excludeRegions = null
if(opts.dgv) {
    DGV dgv = new DGV(opts.dgv).parse() 
    excludeRegions = dgv.grep { Region cnv ->
        
        // Sample sizes that are too small will not contain representative population frequencies
        if(cnv.sampleSize < 5)
            return false
            
        // We are only interested in high frequency CNVs. If we try to exclude all the low frequency ones
        // then we will lose very much (perhaps all) of our genomic regions to simulate with
        float cnvCount = cnv.observedGains + cnv.observedLosses
        float freq = cnvCount / cnv.sampleSize.toFloat()
        if(freq < 0.01) 
            return false
            
        return true
    } as Regions
    println "Regions excluded from simulation are $excludeRegions"
    
    // Debug
    excludeRegions.save("excluded.bed")
}

    
int numRanges = opts.n.toInteger()

if(opts.cov && !opts.abed) {
    System.err.println "Option -cov requires the BED file of covered regions"
    System.exit(1)
}

// Index the bam files we have by sample name
//Map<String,String> bamFiles = samples.grep { entry -> entry.value.files.bam != null }.collectEntries { entry-> [entry.key, entry.value.info.files.bam[0] ] }
Map<String,String> bamFiles = samples.grep { entry -> entry.value.files.bam != null }.collectEntries { entry-> [entry.key, entry.value.files.bam[0] ] }

// For each bam file provided as a command line argument
if(opts.bams) {
	opts.bams.each { bam ->
		bamFiles[new SAM(bam).samples[0]] = bam
	}
}

def missingBams = (females+males).grep { bamFiles[it] == null }
if(missingBams) {
	System.err.println "The following samples specified did not have BAM files available either from meta data or as arguments:\n"
	System.err.println(missingBams.join('\n'))
	System.err.println "\nPlease specify BAM files for these samples using the -bam option"
	System.exit(1)
}

Random random = opts.seed ? new Random(opts.seed.toLong()) : new Random()
cnvs = []
for(female in females) {
    println female.center(100,"=")
    // Parse out the male    
    println "Males are $males"
            
    int selectedMale = Math.floor(random.nextDouble() * males.size())
    
    def selectedMaleSample = males[selectedMale]
    println "Selected male $selectedMale (${selectedMaleSample})"
        
    // Select a region - this is rather simple
    BED bed = new BED(opts.r).load()
        
    CNVSimulator simulator = new CNVSimulator(bamFiles[female], bamFiles[males[selectedMale]])
    simulator.random = random
    if(opts.mode)
        simulator.simulationMode = opts.mode
    
    Region r = simulator.selectRegion(bed, numRanges, excludeRegions)
    println "Seed range for $selectedMaleSample is $r"
    
    String outputFileName = female + "_" + males[selectedMale] + "_$r.from-${r.to}.bam"
    if(opts.d) {
        outputFileName = new File(new File(opts.d),outputFileName).canonicalPath
    }
    
    if(opts.cov) {
        Regions coveredRegions = new BED(opts.abed).load().reduce()
        simulator.setTargetCoverage(coveredRegions, opts.cov.toDouble())
    }
        
	if(opts.threads)
		simulator.concurrency = opts.threads.toInteger()
    
    simulator.createBam(outputFileName, "chrX:$r.from-$r.to")
	
    cnvs.add([female:female,male:males[selectedMale],chr:"chrX",start:r.from, end:r.to, file:outputFileName])
}
    
new File(opts.o).withWriter { w -> 
    w.println cnvs.collect { [it.chr, it.start, it.end, it.file ].join("\t") }.join("\n")
}

if(opts.so) {
    new File(opts.so).withWriter { w ->
        w.println cnvs.collect { cnv ->
            println "Bams from $cnv.female are ${samples[cnv.female].files.bam}"
            samples[cnv.female].files.bam = [cnv.file]
            println "After setting, bams from $cnv.female are ${samples[cnv.female].files.bam}"
            new SampleInfo(sample: cnv.female, batch: new File(opts.d).name, target:"SIM", files: samples[cnv.female].files, sex: "FEMALE").toTsv() 
        }.join("\n")
    }
}
    
println "Wrote CNVs to output file ${opts.o}"
if(opts.so)
    println "Wrote Sample information to output file ${opts.so}"
