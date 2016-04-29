import groovy.lang.Closure;

import java.util.Map;


class ConiferResults extends CNVResults {
    

    ConiferResults(String fileName) {
        super(fileName, 1, 2, 3)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        
        return super.load(options) { Region r ->
            r.sample = r.sampleID
            r.size = r.to - r.from
            r.quality = Math.log10(r.size == 0 ? 1 : r.size)
            if(r.state == "del")
				r.type = 'DEL'
            else
            if(r.state == "dup")
				r.type = 'DUP'
            else
                r.type = 'UNKNOWN'
        }
    }
    
    static void main(String [] args) {
       def con = new ConiferResults("testdata/test.conifer.cnvs.tsv").load()
       
       println con[0].sample + " region " + con[0].toString() + " has " + con[0].type + ' with quality ' + con[0].quality
    }
}
