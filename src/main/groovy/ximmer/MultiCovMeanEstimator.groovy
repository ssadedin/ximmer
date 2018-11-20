package ximmer

import java.io.Reader
import java.util.List
import java.util.Map

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Log

@CompileStatic
@Log
class MultiCovMeanEstimator implements MeanEstimator {
    
    Reader multicovFile

    public MultiCovMeanEstimator(Reader multicovFile) {
        this.multicovFile = multicovFile;
    }

    @Override
    public Map<String, Double> calculateMeans(List<String> samples) {
        // The file is small, so I am just going to read it all and then
        // strip off the lines that are marked NOJSON
        String json = multicovFile.readLines().grep { String line -> !line.endsWith('// NOJSON') }.join('\n')
        
        Map covs = (Map<String,Map>)new JsonSlurper().parseText(json)
        
        Map<String,Double> means = covs.means
        
        log.info "Means found in json cov are: $means"
        
        return samples.collectEntries { [it, means[it] ] };
    }

}
