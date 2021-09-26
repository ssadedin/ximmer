package ximmer.results

class SavvyCNVResults extends CNVResults {
    
    String sourceFile
    String sample
    
    SavvyCNVResults(String sourceFile, String sample=null) {
        super(sourceFile, 0, 1, 2);
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		this.load()
    }
    
    SavvyCNVResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t", columnNames:['chromosome','start','end','type','block_support','block_span','qual','qual_rel','sample']]) { r ->
			r.sample = this.sample.replaceAll('\\.coverageBinner$','')
            r.type = r.type
			r.start = r.start
			r.end = r.end
            r.qual = r.qual
        }
    }
}
