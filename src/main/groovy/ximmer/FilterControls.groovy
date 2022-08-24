package ximmer

import gngs.Cli
import gngs.Utils
import graxxia.Matrix
import graxxia.Stats
import groovy.json.JsonSlurper
import groovy.util.logging.Log

/**
 * A tool that analyses correlation output (see MultiCov) to identify controls where all samples
 * have correlation is smaller than a set threshold.
 * 
 * @author Simon Sadedin
 */
@Log
class FilterControls {

    /**
     * Map keyed on sample id to value of a Map with every sample's correlation
     */
    Map<String,Map<String,Double>> correlations = null

    List<String> controls
    List<String> testSamples

    double threshold
    
    double splitGroupThreshold = 0.82

    double targetCorrelationThreshold = 0.87

    int minPartitionSize = 20

    /**
     * Output directory to write control groups into
     */
    File outputDirectory

    FilterControls(String correlationFile, List<String> controls, double threshold) {
        this.correlations =  new JsonSlurper().parseText(new File(correlationFile).readLines()[1..-1].join('\n'))
        this.controls = controls
        this.threshold = threshold
        this.testSamples = this.correlations*.key - controls
    }

    void run() {

        log.info "Total sample set is ${correlations*.key}"
        log.info "Filtering ${controls.size()} controls to correlation ${threshold}"

        log.info "${testSamples.size()} non-control samples are: $testSamples"
        
        List<String> filteredControls = controls.grep {
            isAcceptableControl(it)
        } 
        
        List<SampleSet> controlGroups = createControlGroups(filteredControls)
        
        log.info "Selected controls are: ${filteredControls}"
        
        println testSamples.join('\n') + '\n' + filteredControls.join('\n') + '\n'
        
        controlGroups.eachWithIndex { group, i ->
            
            Stats stats = group.getCorrelationStats()
            
            File outputFile = new File(outputDirectory, "control_set_${i}.txt")
            outputFile.text = group.samples.keySet().join('\n') + '\n'
            log.info "Wrote ${group.samples.size()} samples for control set ${i} with minCorr=${String.format('%.2f',stats.min)}, meanCorr=${String.format('%.2f',stats.mean)} to ${outputFile}"
        }
    }
    
    List<SampleSet> createControlGroups(List<String> filteredControls) {
        
        List<String> allSamples = testSamples + filteredControls
        
       
        SampleSet initialSet = new SampleSet(samples:correlations.collectEntries { sample, correlations ->
            [ sample, new Sample(sample:sample, correlations:correlations)]
        })
        
        List results = [initialSet]
        
        Sample worstSample = initialSet.worst()
        if(worstSample.worst.value > splitGroupThreshold) {
            log.info "All samples are well matched, no requirement to split control sets"
        }
        else {
            log.info "Worst sample correlation for $worstSample ($worstSample.worst) is below threshold $splitGroupThreshold"
            results = initialSet.splitSet(worstSample.sample, targetCorrelationThreshold)
            
            log.info "Partitioned to raw sets of size ${results*.samples*.size()}"
            
            // Add back controls that have good correlation
            int partitionIndex = 0
            for(SampleSet partition in results) {
                ++partitionIndex
                filteredControls.each { String controlId ->
                    Sample control = initialSet.samples[controlId]
                    if(!partition.samples.containsKey(controlId) && partition.samples.every { correlations[it.key][controlId] > targetCorrelationThreshold }) {
                        
                        log.info "Control $controlId can be rescued into partition $partitionIndex"

                        partition.addSample(control)
                    }
                }
                
                // What if we still do not have enough good controls in the partition?
                // Then top up with the best we can from the remainder of the control set
                if(partition.samples.size() < minPartitionSize)  {
                    
                    log.info "Partition $partitionIndex is still too small, adding in supplementary control samples"

                    // Find the best remaining controls and add it in
                    List bestControls = findBestSuboptimalControls(partition, worstSample)
                    
                    log.info "The best controls to add to partition $partitionIndex are: $bestControls"
                    bestControls.each { controlId ->
                        partition.addSample(initialSet.samples[controlId])
                    }
                }
            }
        }

        return results
    }

    /**
     * Search for controls not already included in the target sample set that are best matched
     * to the given key sample.
     * 
     * @param targetSampleSet   sample set into which controls are to be inserted
     * @param keySample         sample to which to attempt to match the controls
     * @return
     */
    private List findBestSuboptimalControls(SampleSet targetSampleSet, Sample keySample) {
        // Sort controls by their correlation to worstSample
        List bestControls =
                controls
                .grep { !targetSampleSet.samples.containsKey(it) }
                .sort { controlId ->
                    -correlations[controlId][keySample.sample]
                }
                .take(minPartitionSize - targetSampleSet.samples.size())
        return bestControls
    }
    
    
    /**
     * Test if a control has better than the given threshold correlation to at least one test sample
     * 
     * @param control
     * @return
     */
    boolean isAcceptableControl(String control) {
        int numTestSamples = testSamples.count {  String testSample ->
            this.correlations[testSample][control] > this.threshold
        }
        
        if(numTestSamples>0) {
            log.info "Control $control is well matched to ${numTestSamples} test samples"
            return true
        }
        else {
            log.info "Control $control is not well matched to any test samples: this control will be rejected"
            return false
        }
    }

    static void main(String [] args) {

        Utils.configureSimpleLogging()

        Cli cli = new Cli(usage: 'FilterControls <options>')
        cli.with {
            corr 'JSON file containing samplewise correlations', args:1, required: true
            control 'Specify a control sample', args: Cli.UNLIMITED, required: true
            thresh 'Minimum correlation threshold to accept a control (0.9)', args: 1, required: false
            splitThreshold 'correlation threshold at which to split test sample to a separate group', args:1, required: false, type: Double
            minimumGroupSize 'Minimum target number of samples to achieve in each group. If a group contains less samples then poor-matching controls may be added to increase the contol set size', required: false, args: 1, type: Integer
            o 'Output directory to write multiple control sets into', longOpt: 'outputDirectory', args:1, required: false
        }

        OptionAccessor opts = cli.parse(args)
        if(!opts)
            System.exit(1)

        double threshold = opts.thresh? opts.thresh.toDouble() : 0.9
        def fc = new FilterControls(opts.corr, opts.controls, threshold)
        
        if(opts.o) {
            fc.outputDirectory = new File(opts.o)
        }
        
        if(opts.splitThreshold)
            fc.splitGroupThreshold = opts.splitThreshold
        
        if(opts.minimumGroupSize)
            fc.minPartitionSize = opts.minimumGroupSize

        fc.run()
    }
}

class Sample {
    String sample
    
    Map<String,Double> correlations
    
    Map.Entry<String,Double> getWorst() {
        correlations.min { it.value }
    }
    
    String toString() { "Sample $sample"}
}

class SampleSet {
   Map<String, Sample> samples
   
   /**
    * Add the given other sample to this sample set.
    * 
    * Adjusts the correlation maps of the existing samples to include the control
    * and adds the control itself as a new {@link Sample} object with correlations
    * narrowed to the existing samples in this set.
    * 
    * @param control
    */
   void addSample(Sample sampleToAdd) {
       
       String controlId = sampleToAdd.sample
       
        def controlCorrelations = this.samples.keySet().collectEntries { s ->
           [s, sampleToAdd.correlations[s]]
        }
        controlCorrelations[controlId] = 1.0d
        
        this.samples.each { String id, Sample sample ->
            sample.correlations[controlId] = sampleToAdd.correlations[id]
        }

        this.samples[controlId] = 
            new Sample(sample:controlId, correlations: controlCorrelations)
   }
       
    
    Stats getCorrelationStats() {
        Stats.from(samples.collectMany { it.value.correlations*.value })
    }
    
    Sample worst() {
        return samples.min { e ->
            e.value.correlations*.value.min()
        }.value
    }
    
    List<SampleSet> splitSet(String seed, Double targetCorrelation) {
        
        List<Map<String,Sample>> partitions = samples.split { e ->
            e.value.correlations[seed] > targetCorrelation
        }*.collectEntries()
        .collect { Map<String,Sample> partitionSamples ->
            partitionSamples.collectEntries { e ->
                Sample sample = e.value
                def matchingCorrelations = sample.correlations.grep { it.key in partitionSamples}.collectEntries()
                [sample.sample, new Sample(sample: sample.sample, correlations: matchingCorrelations)]
            }
        }
       
        return partitions.collect { p ->
            new SampleSet(samples: p)
        }
    }
    
    String toString() { "${samples.size()} samples min corr=${String.format('%.2f',this.worst().worst.value)}: ${samples*.key.join(', ')}" }
}

