package ximmer.results

class ParallaxResults extends CNVResults {
    
    String sourceFile
    String sample
    
    public ParallaxResults(String sourceFile) {
        super(sourceFile, 0, 1, 2);
        this.separator = "\t"
        this.sourceFile = sourceFile
        this.sample = sourceFile.replaceAll('\\..*$','')
        List parts = sample.tokenize('_')
        
        // VCGS specific logic
        // TODO: fix to make sample passable as param
        if(parts == 6) {
            this.sample = parts[4]
        }
    }
    
    ParallaxResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t", columnNames:['chr','start','end','quality']]) { r ->
            r.type = "DEL"
        }
    }
}
