package ximmer.results

class ParallaxResults extends CNVResults {
    
    String sourceFile
    String sample
    
    ParallaxResults(String sourceFile, String sample=null) {
        super(sourceFile, 0, 1, 2);
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		this.load()
    }
    
    ParallaxResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t", columnNames:['chr','start','end','type','quality']]) { r ->
			r.sample = this.sample
            r.type = r.type
			r.start = r.from
			r.end = r.to
        }
    }
}
