package ximmer.results

class SavvyCNVResults extends CNVResults {
    
    String sourceFile
    String sample

    final static CN_TO_TYPE = [ 
            "Deletion": "DEL", 
            "Duplication" :"DUP",
    ]
 
    
    SavvyCNVResults(String sourceFile, String sample=null) {
        super(sourceFile, 0, 1, 2);
        this.sourceFile = sourceFile
        this.sample = sample ? sample : getSampleFromFile(sourceFile)
        this.load()
    }
    
    SavvyCNVResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t"]) { r ->
            r.sample = new File(r.sample).name.replaceAll('\\.coverageBinner$','')
            r.type = CN_TO_TYPE[r.type] ?: 'OTHER'
            r.quality = r.qual
        }
    }
}
