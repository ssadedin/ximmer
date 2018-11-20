package ximmer

import java.util.Map
import java.util.List

import gngs.*

/**
 * A relatively expensive but accurate estimator of mean coverage
 * which calculates coverage statistics by actually computing the read 
 * pileup across the target regions0
 * 
 * @author Simon Sadedin
 */
class PileupMeanEstimator implements MeanEstimator {
    
    Regions meanRegions
    
    Map<String,SAM> bams

    public PileupMeanEstimator(Regions meanRegions, Map<String, SAM> bams) {
        super();
        this.meanRegions = meanRegions;
        this.bams = bams;
    }

    @Override
    public Map<String, Double> calculateMeans(List<String> samples) {
        Map<String, Double> result  = Collections.synchronizedMap([samples, 
                                     samples.collectParallel { bams[it].coverageStatistics(meanRegions).mean }]
                                               .transpose().collectEntries()) 
       return result
    }

}
