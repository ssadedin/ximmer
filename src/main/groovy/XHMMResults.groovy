
import com.xlson.groovycsv.PropertyMapper

import gngs.RangedData
import gngs.Region


class XHMMResults extends CNVResults {
    
    XHMMResults(String fileName) {
        super(fileName, -1, -1, -1)
    }
    
    @Override
    public RangedData load(Map options, Closure c=null) {
        return super.load(options) { Region r ->
            r.sample = r.SAMPLE
            r.type = r.CNV
            r.quality = r.Q_SOME
        }
    } 
    
	@Override
    protected Region parseRegion(PropertyMapper line) {
        return new Region(line[2])
    }


    static void main(String [] args) {
        def xhmm = new XHMMResults("testdata/asdsibs.params.xhmm_discover.xcnv").load()
        println xhmm[0].toString() + " sample " + xhmm[0].sample
       
    }
}
