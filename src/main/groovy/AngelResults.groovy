
import java.util.Map;


class AngelResults extends CNVResults {
    
    AngelResults(String fileName) {
        super(fileName)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        return super.load(options + [columnNames: ['chr','start','end','sample','quality']]) { r ->
            r.type = "DEL"
        }
    }
    
    static void main(String [] args) {
       def angel = new AngelResults("testdata/autosomes.angelhmm.cnvs.bed").load()
       
       println angel[0].sample + " region " + angel[0].toString() + " has " + angel[0].type + " and quality " + angel[0].quality
    }
}
