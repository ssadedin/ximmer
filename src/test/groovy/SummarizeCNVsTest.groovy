import static org.junit.Assert.*

import org.junit.Test
import org.junit.*
import gngs.*

import ximmer.results.*


class SummarizeCNVsTest {
	
	Map<String,RangedData> results = [:]
	
    SummarizeCNVs scnvs = new SummarizeCNVs()
    
    Region cnv = new Region('chr1:1000-2000')
    
	@Test
	public void 'basic CNV annotation'() {
        cnv.sample = 'FOO'
        
        Map data = scnvs.cnvToMap([], [], SummarizeCNVs.DEFAULT_JS_COLUMNS + [], cnv)
           
        Utils.table([data])
        
        assert data.sample == 'FOO'
	}
    
    @Test
    void 'test caller span annotations'() {
        List callers = ['ed','xhmm']
        
        Region ed = new Region('chr1:900-1300')
        ed.quality = 10
        ed.calls = [ed]
        cnv.ed = ed
        
        Map data = scnvs.cnvToMap(callers, [], scnvs.computeColumns(callers,[]), cnv)
        
        println data
        
        Utils.table([data])
        
        assert data.ed == 'TRUE'
        assert data.xhmm == 'FALSE'
        assert data.calls.ed[0] == [900,1300,10]
    }
}
