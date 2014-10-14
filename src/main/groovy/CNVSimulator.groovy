import groovy.transform.CompileStatic;
import groovyx.gpars.GParsPool;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTagUtil;

/**
 * Produce a BAM file that simulates the presence of a heterozygous deletion in
 * the X chromosome of a given female sample, by substituting reads from a male 
 * sample in a specified region.
 * <p>
 * A number of issues need to be dealt with to make this realistic:
 * 
 * <li>The coverage of the two samples must be normalised to account for the library size
 * <li>Breaking the coverage in a region where there are reads will cause a discontinuity.
 *     Therefore it is important to detect the regions of no coverage so as to 
 *     perform the "swap" in a location where no artefactual read discontinuities will occur.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class CNVSimulator {
    
    /**
     * Chromosomes that will be ignored for calculates that are only performed
     * on diploid chromosomes
     */
    public static List<String> NON_AUTOSOMES = ["chrX",  "chrY", "chrM"]
    
    SAM maleBam 
    
    SAM femaleBam 
    
    Regions maleReadRegions 
    
    Regions femaleReadRegions
    
    /**
     * Fraction of female reads to use from X chromosome
     * Value < 0 means not calculated yet.
     */
    double femaleDownSampleRate = -1.0d
    
    /**
     * Fraction of male reads to use from region to be imported to female 
     * alignment
     * Value < 0 means not calculated yet.
     */
    double maleDownSampleRate = -1.0d
    
    /**
     * Maximum number of attempts to select a random region satisfying all constraints
     */
    int maxSelectionAttempts = 20
    
    /**
     * Default concurrency used for some multithreadded parts
     */
    int concurrency = 2
    
    /**
     * The kind of simulation to do. One of "replace" or "downsample"
     */
    String simulationMode="replace"
    
    /**
     * The mean coverage to be targeted in the output samples.
     * This value is an optional input parameter. If provided,
     * the coverage.
     */
    private double targetCoverage = -1.0d
    
    /**
     * If targetCoverage is specified, the regions over which it is calculated
     */
    private Regions targetRegions = null
    
    /**
     * Random number generator to use for downsampling. The user can set it to
     * predetermine the seed or use different sampling strategy.
     */
    Random random = null
    
    public CNVSimulator(String femaleBam, String maleBam) {
        this.maleBam = new SAM(maleBam)
        this.femaleBam = new SAM(femaleBam)
        this.random = new Random()
    }
    
    void setTargetCoverage(Regions targetRegions, double value) {
        this.targetRegions = targetRegions
        this.targetCoverage = value
        this.random = new Random()
    }
    
    static void main(String [] args) {
        
        println "=" * 100
        println "CNV Simulator " + new Date().toString()
        println "=" * 100
        
        Cli cli = new Cli(usage:"CNVSimulator [options]")
        cli.with {
            f "female sample", longOpt: "female", args:1, required:true
            m "male sample", longOpt: "male", args:1, required:true
            t "Number of threads to use", args:1 
            region "region to turn into a CNV", args:1, required:true
            cov "target coverage for output samples (optional)", args:1
            bed "BED file describing covered regions (required for -cov), eg: Amplicon bed file for HaloPlex", args:1
            o "output file", args:1, required:true
            seed "Random seed, specify to make simulations reproducible", args:1
            mode "Method for simulating deletions: downsampling reads=downsample, replacement from male=replace, both=use both methods (2 output files)", args:1
        }
        
        def opts = cli.parse(args)
        if(!opts) 
            System.exit(1)
            
        def simulator = new CNVSimulator(opts.f, opts.m)
        
        if(opts.t)
            simulator.concurrency = opts.t.toInteger()
            
        if(opts.seed) 
            simulator.random = new Random(opts.seed.toLong())
        
        if(opts.cov) {
            if(!opts.bed) {
                System.err.println "Please provide option -bed with the amplicon BED file, to specify regions to calculate coverage over"
                System.exit(1)
            }
            simulator.setTargetCoverage(new BED(opts.bed).load(), opts.cov.toDouble())  
        }
        
        def simulationModes = ["replace","downsample"]
        if(opts.mode) { 
            if(!(opts.mode in simulationModes)) {
                System.err.println "Please specify one of " + simulationModes + " for 'mode' option"
                System.exit()
            }
                
            simulator.simulationMode = opts.mode
        }
        
        simulator.createBam(opts.o, opts.region)
    }
    
    @CompileStatic
    int countAutosomeReads(SAM bam) {
        int count=0
        bam.eachRecord { SAMRecord r ->
            if(NON_AUTOSOMES.contains(r.referenceName))
                return
            ++count    
        }
        return count
    }
    
    void calculateDownsampleRates() {
        
        // Default is not to downsample - which one we downsample
        // depends on which one is higher to start with
        femaleDownSampleRate = 1.0d
        maleDownSampleRate = 1.0d
        
        GParsPool.withPool(this.concurrency) {
            println "Using $concurrency threads"
            if(this.targetCoverage > 0.0d) {
		        Regions autosomalRegions = null
		        synchronized(this) { // cause write flush after we build this block
		            autosomalRegions = this.targetRegions.reduce().grep { !NON_AUTOSOMES.contains(it.chr) }
		        }
	        
                List femaleStats = null
                List maleStats = null
                println "Determining coverage of female autosomal regions ..."
                ProgressCounter counter = new ProgressCounter()
                femaleStats = autosomalRegions.collectParallel { counter.count(); femaleBam.coverageStatistics(it.chr, it.from, it.to); }
                
                println "Determining coverage of male autosomal regions ..."
                maleStats = autosomalRegions.collectParallel { counter.count(); maleBam.coverageStatistics(it.chr, it.from, it.to); }
                
				double femaleSum = (femaleStats.sum { it.mean * it.n })
                double femaleMean = femaleSum / femaleStats.sum { it.n }
                double maleMean = maleStats.sum { it.mean * it.n } / maleStats.sum { it.n }
                
                this.femaleDownSampleRate = this.targetCoverage / femaleMean
                this.maleDownSampleRate = this.targetCoverage / maleMean
                
                println "Female coverage = $femaleMean: downsampling reads by $femaleDownSampleRate"
                println "Male coverage = $maleMean: downsampling reads by $maleDownSampleRate"
            }
            else {
                
                int maleReadCount = countAutosomeReads(maleBam)
                int femaleReadCount = countAutosomeReads(femaleBam)
                
                
                if(maleReadCount < femaleReadCount) {
                    // Female sample is too high, so downsample the female sample
                    femaleDownSampleRate = ((float)maleReadCount) / ((float)femaleReadCount)
                }
                else {
                    // Male sample is too high, so downsample the male sample
                    maleDownSampleRate = ((float)femaleReadCount) / ((float)maleReadCount)
                }
                
                println "Male reads ($maleBam.samFile)  " + maleReadCount + " downsample rate = " + maleDownSampleRate
                println "Female reads ($femaleBam.samFile)  " + femaleReadCount + " downsample rate = " + femaleDownSampleRate
            }
        }
    }
    
    void createBam(String outputFileName, String region) {
        
        calculateDownsampleRates()
        
        if(this.simulationMode=="replace") {
            this.createBamByReplacement(outputFileName, region)
        }
        else
        if(this.simulationMode=="downsample") {
            this.createBamByDownSample(outputFileName, region)
        }
        else
        if(this.simulationMode=="both") {
            this.createBamByReplacement(outputFileName, region)
            
            def bothName = outputFileName.replaceAll('.bam$','.downsample.bam')
            this.createBamByDownSample(bothName, region)
            println "Created $bothName"
        }        
        else {
            throw new IllegalArgumentException("Simulation mode of $simulationMode is not recognized. Please use one of [replace,downsample,both]")
        }
        
        println "Created $outputFileName"
    }
    
    void createBamByReplacement(String outputFileName, String region) {
        
        println "Creating $outputFileName using replace mode"
        
		Region cleanRegion = new Region(region)
        
        println "Querying source reads from Male alignment ..."
        List maleRegionReads = []
        int count = 0
        maleBam.eachPair(cleanRegion.chr,cleanRegion.from,cleanRegion.to) { r1, r2 ->
            maleRegionReads << [ r1, r2]
            ++count
        }
        
        println "Found $count read pairs in region over male alignment " + this.maleBam.samples[0]
        if(count == 0) 
            throw new RuntimeException("Male alignment has no reads over targeted region for CNV")
            
        if(this.random == null)
            random = new Random()
        
        final String chr = cleanRegion.chr
		final int start = cleanRegion.from
		final int end = cleanRegion.to
        
        final String rgId = femaleBam.samFileReader.fileHeader.getReadGroups()[0].getId()
        
        // Read each BAM 
        femaleBam.withWriter(outputFileName) { SAMFileWriter  writer ->
            femaleBam.eachPair { SAMRecord r1, SAMRecord r2 ->
                
                if(r1 == null|| r2 == null)
                    return
                
			    def boundaries = [r1.alignmentStart, r1.alignmentEnd, r2.alignmentStart, r2.alignmentEnd]
				int rStart = boundaries.min()
				int rEnd = boundaries.max()
				
                if(cleanRegion.overlaps(chr, rStart, rEnd))
                    return
                
                // For regions outside that which we are simulating for,
                // downsample the reads to make the mean coverage
                // the same as that from the lower coverage file
                // (the lower coverage file might be this one or the other file that is the source for simulated reads)
                if(random.nextFloat() < femaleDownSampleRate) {
                    writer.addAlignment(r1)
                    writer.addAlignment(r2)
                }
            }
            
            for(List<SAMRecord> r in maleRegionReads) {
                if(random.nextFloat() < maleDownSampleRate) {
					r[0].setAttribute("DL", "1")
                    r[0].setAttribute(SAMTagUtil.getSingleton().RG, rgId);
                    writer.addAlignment(r[0])
					r[1].setAttribute("DL","1")
                    r[1].setAttribute(SAMTagUtil.getSingleton().RG, rgId);
                    writer.addAlignment(r[1])
                }
            }
        }
    }
    
    /**
     * This method creates a BAM file with simulated deletions purely by randomly selecting
     * only 50% of the reads over the CNV region.
     * 
     * @param outputFileName
     * @param region
     */
    void createBamByDownSample(String outputFileName, String region) {
        
        println "Creating $outputFileName using downsample mode"
        
        Region cleanRegion = new Region(region)
        
        // Inside regions of simulated deletions, downsample to 50% of the otherwise estimated rate
        // technical note: the femaleDownSampleRate takes into account the coverage in the supplied
        //                 male alignment. This is necessary to make the output comparable
        //                 in terms of coverage when using 'both'.
        double deletionDownSampleRate = 0.5d * femaleDownSampleRate
        
        String chr = cleanRegion.chr
        
        // Read each BAM 
        femaleBam.withWriter(outputFileName) { SAMFileWriter  writer ->
            femaleBam.eachPair { SAMRecord r1, SAMRecord r2 ->
                
                if(r1 == null|| r2 == null)
                    return
                
			    def boundaries = [r1.alignmentStart, r1.alignmentEnd, r2.alignmentStart, r2.alignmentEnd]
				int rStart = boundaries.min()
				int rEnd = boundaries.max()
                
                float downSampleRate = femaleDownSampleRate
				
                if(cleanRegion.overlaps(chr, rStart, rEnd)) {
                    downSampleRate = deletionDownSampleRate
                }
                
                // For regions outside that which we are simulating for,
                // downsample the reads to make the mean coverage
                // the same as that from the lower coverage file
                // strictly not necesary for a pure downsample simulation, but we want to make
                // output that is comparable with that produced when reads are sourced from a male
                if(random.nextFloat() < downSampleRate) {
                    writer.addAlignment(r1)
                    writer.addAlignment(r2)
                }
            }
        }        
    }
	
	/**
	 * Expand the given region until it is "clean" on either side in both the
	 * male and female files.
	 * 
	 * @param region
	 * @return
	 */
	Region findCleanRegion(Region seedRegion) {
		
        if(maleReadRegions == null)
            maleReadRegions = maleBam.toPairRegions("chrX",0,0,500)
        
		println "Overlaps of male region with target are " + maleReadRegions.getOverlaps(seedRegion).collect { it.from + "-" + it.to }
		
        if(femaleReadRegions == null)
            femaleReadRegions = femaleBam.toPairRegions("chrX",0,0,500)
            
		println "Overlaps of female region with target are " + femaleReadRegions.getOverlaps(seedRegion).collect { it.from + "-" + it.to }
        
		Regions combinedRegions = maleReadRegions.reduce()
		femaleReadRegions.reduce().each { r ->
			combinedRegions.addRegion(r)
		}
		
		combinedRegions = combinedRegions.reduce()
		
		// Finally, find the overlaps with the desired region
		List<IntRange> overlaps = combinedRegions.getOverlaps(seedRegion)
	
		// We will start from the end points of this range
		Region result = new Region("chrX", overlaps*.from.min()..overlaps*.to.max())
        
        println "Original region $seedRegion expanded to $result [female=" + this.femaleBam.samples[0] + ", male=" + this.maleBam.samples[0] + "]"
		
		return result
	}
    
    /**
     * Select a candidate region to simulate a deletion from a set of supplied regions, 
     * while avoiding a given set of excluded regions.
     * 
     * @param fromRegions   target regions to consider as candidates for simulating a deletion
     * @param numRanges     number of contiguous ranges to include from fromRegions
     * @param excludRegions regions to exclude from consideration as targets
     * @return
     */
    Region selectRegion(Regions fromRegions, int numRanges, Regions excludeRegions=null) {
        Region cleanRegion  = null
        int attemptCount = 0
        while(true) {
            int selectedRange =  (int)Math.floor(random.nextDouble() * (fromRegions.allRanges["chrX"].size()-numRanges)) 
            List<Range> regions = fromRegions.allRanges["chrX"][selectedRange..(selectedRange+numRanges-1)] 
            Range r = (regions[0].from)..(regions[-1].to)
            
            Region seedRegion = new Region('chrX',r)
            
            // The seed region is directly derived from the target region
            // The problem is, many capture technologies will capture reads to either side of the
            // target region. This means we need to expand the region either side until we observe 
            // no reads, so that a "clean" swap can be made of reads without causing any artefacts.
            cleanRegion = this.findCleanRegion(seedRegion)
            
            // Check if the cleaned up region overlaps any excluded regions
            if(!excludeRegions || !excludeRegions.overlaps(cleanRegion))
                break
               
            ++attemptCount
            if(attemptCount > 20)
                throw new RuntimeException("Failed to identify a non-excluded region for placement of deletion after ")
                
            println "Selected range $r.from-$r.to overlaps one or more excluded regions: trying again"
        }
        return cleanRegion
    }
}
