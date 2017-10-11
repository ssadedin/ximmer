import groovy.lang.Closure;

import java.util.Map;


class CodexResults extends CNVResults {
    

    CodexResults(String fileName) {
        super(fileName, 1, 3, 4)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        
        return super.load(options) { Region r ->
            r.sample = r.sample_name
            r.size = r.to - r.from
            r.quality = r.lratio
            if(r.cnv == "dup")
				r.type = 'dup'
            else
            if(r.cnv == "del")
				r.type = 'DEL'
            else
                r.type = 'UNKNOWN'
        }
    }
    
    static void main(String [] args) {
       def con = new CodexResults("codex.tsv").load()
       
       println con[0].sample + " region " + con[0].toString() + " has " + con[0].type + ' with quality ' + con[0].quality
    }
}
