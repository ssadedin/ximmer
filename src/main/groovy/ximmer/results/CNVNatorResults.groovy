package ximmer.results

import gngs.*

class CNVNatorResults extends CNVResults {
    
    String sourceFile
    String sample
	String type
    
    public CNVNatorResults(String sourceFile, String sample=null) {
		super(sourceFile)
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		
		def v = VCF.parse(this.sourceFile) {
            
            if(it.filter != 'PASS')
                return false
            
			Region r = new Region(it.chr, it.pos, it.info.END.toInteger())
			r.type = it.info.SVTYPE
			r.sample = this.sample
			r.start = it.pos
			r.end = it.info.END.toInteger()
			r.quality = it.info['natorQ0'].toDouble()
			addRegion(r)
            
            return false
		}
    }
	
}

