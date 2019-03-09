package ximmer.results

import gngs.*

class LumpyResults extends CNVResults {
    
    String sourceFile
    String sample
    
    public LumpyResults(String sourceFile, String sample=null) {
		super(sourceFile)
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		
		VCF.parse(this.sourceFile) { v ->
            if(!(v.info.SVTYPE in ['DEL','DUP'])) {
                return false
            }
            
            if(Region.isMinorContig(v.chr))
                return false
                
			Region r = new Region(v.chr, v.pos, v.info.END.toInteger())
			r.type = v.info.SVTYPE
			r.sample = this.sample
			r.start = v.pos
			r.end = v.info.END.toInteger()
			r.quality = [v.info.PE?:0,v.info.SU?:0, v.info.SR?:0]*.toInteger().sum()
            
			addRegion(r)
            
            return false
		}
    }
	
}

