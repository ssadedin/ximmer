package ximmer.results

import gngs.*

class DellyResults extends CNVResults {
    
    String sourceFile
    String sample
    
    public DellyResults(String sourceFile, String sample=null) {
		super(sourceFile)
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		
		VCF.parse(this.sourceFile) {
            if(!(it.info.SVTYPE in ['DEL','DUP'])) {
                return
            }
                
			Region r = new Region(it.chr, it.pos, it.info.END.toInteger())
			r.type = it.info.SVTYPE
			r.sample = this.sample
			r.start = it.pos
			r.end = it.info.END.toInteger()
			r.quality = it.filter == 'PASS' ? 20 : 0
            
			addRegion(r)
            
            return false
		}
    }
	
}

