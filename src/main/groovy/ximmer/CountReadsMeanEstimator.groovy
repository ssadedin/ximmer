package ximmer

import java.util.List
import java.util.Map

import gngs.*
import groovy.util.logging.Log

@Log
class CountReadsMeanEstimator implements MeanEstimator {
    
    Regions targetRegions
    
    Map<String,SAM> bams
    
    public CountReadsMeanEstimator(Regions targetRegions, Map<String, SAM> bams) {
        this.targetRegions = targetRegions;
        this.bams = bams;
    }

    @Override
    public Map<String, Double> calculateMeans(List<String> samples) {
        
        int targetSize = targetRegions.size()
        Map<String, Double> means = Collections.synchronizedMap([samples, samples.collectParallel {
                log.info "Counting reads for $it";
                if(!bams[it])
                    System.err.println "No BAM file provided for $it"

                long count=0; bams[it].eachRecord { count += it.readLength };
                count/(double)targetSize
            }
        ].transpose().collectEntries())
    }

}
