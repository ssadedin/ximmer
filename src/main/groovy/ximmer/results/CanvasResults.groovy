package ximmer.results

import gngs.*

class CanvasResults extends CNVResults {
    
    String sourceFile
    String sample
    
    static final Map<String,String> CANVAS_TYPE_MAP = [
        LOSS : 'DEL',
        GAIN : 'DUP'
    ]
    
    public CanvasResults(String sourceFile, String sample=null) {
		super(sourceFile)
        this.sourceFile = sourceFile
		this.sample = sample ? sample : getSampleFromFile(sourceFile)
		
		VCF.parse(this.sourceFile) { Variant v ->
            
            String type = CANVAS_TYPE_MAP[v.id.tokenize(':')[1]]
            if(type == null)
                return false
                
            if(v.quality < 2)
                return
               
			Region r = new Region(v.chr, v.pos, v.info.END.toInteger())
			r.type = type
			r.sample = this.sample
			r.start = v.pos
			r.end = v.info.END.toInteger()
			r.quality = v.qual
            
			addRegion(r)
            return false
		}
    }
	
}

