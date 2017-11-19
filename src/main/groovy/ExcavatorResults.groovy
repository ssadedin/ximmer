import groovy.lang.Closure;

import java.util.Map;

import gngs.Region


class ExcavatorResults extends CNVResults {
    

    ExcavatorResults(String fileName) {
        super(fileName, 1, 2, 3)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        
        return super.load(options + [columnNames: ['sample','chr','start','end', 'ratio','change'], readFirstLine:false]) { Region r ->
            r.size = r.to - r.from
            r.quality = Math.log10(r.size == 0 ? 1 : r.size)
            if(r.type < 0)
				r.type = 'DEL'
            else
            if(r.type > 0)
				r.type = 'DUP'
            else
                r.type = 'UNKNOWN'
        }
    }
    
    static void main(String [] args) {
       def ex = new ExcavatorResults("testdata/example_exacavator.cnvs.tsv").load()
       
       println ex[0].sample + " region " + ex[0].toString() + " has " + ex[0].type + ' with quality ' + ex[0].quality
    }
}
