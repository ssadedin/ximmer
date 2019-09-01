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
                return false
            }
            
            if(it.filter != 'PASS')
                return false
            
            if(Region.isMinorContig(it.chr))
                return false
                
			Region r = new Region(it.chr, it.pos, it.info.END.toInteger())
			r.type = it.info.SVTYPE
			r.sample = this.sample
			r.start = it.pos
			r.end = it.info.END.toInteger()
			r.quality = it.genoTypes[0].GQ.toDouble()
            
			addRegion(r)
            
            return false
		}
    }
	
}

