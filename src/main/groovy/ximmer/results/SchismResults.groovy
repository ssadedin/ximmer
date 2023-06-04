package ximmer.results

import gngs.*
import groovy.json.JsonSlurper

class SchismResults extends CNVResults {
    
    String sourceFile
    List<String> samples
    
    public SchismResults(String sourceFile) {
		super(sourceFile)
        this.sourceFile = sourceFile
        
        
        // Problem: how to deal with the fact that schism puts multiple sample's results
        // into a single file?
        
        // a) assume there is a primary sample for the file and store under that?
        
        def json = new JsonSlurper().parse(new File(sourceFile))
        
        /* Structure of a Schism result:
         * {
                "chr": "chr1",
                "start": 19769913,
                "end": 19769914,
                "sample": {
                    "NA12878": {
                        "obs": 5,
                        "motif": "TTGCACACATGCTCA",
                        "startClips": 5,
                        "endClips": 0
                    }
                },
                "depth": 5,
                "sample_count": 1,
                "cscore": "1",
                "partner": "chr21:10432380",
                "genes": [
                    "TMCO4"
                ],
                "cdsdist": [
                    629
                ],
                "samples": []
            }
         */
        for(Map result in json) {
            def  nonDataKeys = ['depth','sample_count','cscore','partner','cdsdist']
            result.sample.each { sample, data ->
                Region r = new Region(result.chr, result.start, result.end)
                r.type = "BND"
                r.quality = data.obs
                r.sample = sample
                r.details = nonDataKeys.collectEntries { [it, result[it]] } + data
                addRegion(r)
            }
        }
    }
	
}

