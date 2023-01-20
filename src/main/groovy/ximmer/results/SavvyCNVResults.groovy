package ximmer.results

import gngs.Region
import gngs.Regions

class SavvyCNVResults extends CNVResults {
    
    String sourceFile
    String sample
    
    Regions targetRegions

    final static CN_TO_TYPE = [ 
            "Deletion": "DEL", 
            "Duplication" :"DUP",
    ]
 
    
    SavvyCNVResults(String sourceFile, Regions targetRegions, String sample=null) {
        super(sourceFile, 0, 1, 2);
        this.sourceFile = sourceFile
        this.sample = sample ? sample : getSampleFromFile(sourceFile)
        this.targetRegions = targetRegions
        this.load()
    }
    
    SavvyCNVResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t"]) { Region r ->
            
            // Because SavvyCNV does not confine calling to the target regions we want to adjust its
            // calls by intersecting them down to the nearest target region boundary to make them
            // consistent with how other callers are treated
            if(this.targetRegions != null) {
                List<IntRange> overlaps = targetRegions.getOverlaps(r)
                if(overlaps.isEmpty())
                    return false
                    
                int start = overlaps[0].from
                int end = overlaps[-1].to
                r.setRange(start..end)
            }
            
            r.sample = new File(r.sample).name.replaceAll('\\.coverageBinner$','')
            r.type = CN_TO_TYPE[r.type] ?: 'OTHER'
            r.quality = r.qual
        }
    }
}
