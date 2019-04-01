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
        ed.all = [ed]
        ed.best= ed
        cnv.ed = ed
        cnv.xhmm = [:]
        
        Map data = scnvs.cnvToMap(callers, [], scnvs.computeColumns(callers,[]), cnv)
        
        println data
        
        Utils.table([data])
        
        assert data.ed == 'TRUE'
        assert data.xhmm == 'FALSE'
        assert data.calls.ed[0] == [900,1300,10]
    }
    
    @Test
    void 'test single caller are annotated correctly'() {
        Region edCall = new Region('chr1', 1000..2000, quality: 100, sample: 'MrBoo', )
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // Should find 'ed' overlaps
        Region cnv = new Region('chr1', 1200..1800)
        assert scnvs.annotateCaller('MrBoo', cnv,  'ed')
        assert cnv.ed.best.is(edCall)
        assert cnv.ed.supporting.size() == 1
        assert cnv.ed.all.size() == 1
    }
    
    @Test
    void 'test single call only supported when mutal overlap'() {
        Region edCall = new Region('chr1', 1200..1250, quality: 100, sample: 'MrBoo', )
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // Should find 'ed' overlaps
        Region cnv = new Region('chr1', 1000..2000) // tiny overlap
        assert !scnvs.annotateCaller('MrBoo', cnv,  'ed')
        assert cnv.ed.best == null
        assert cnv.ed.supporting.size() == 0
        assert cnv.ed.all.size() == 1
    } 
    
    @Test
    void 'test single call only supported when correct sample'() {
        Region edCall = new Region('chr1', 1200..1250, quality: 100, sample: 'MsFoo', )
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // Should find 'ed' overlaps
        Region cnv = new Region('chr1', 1000..2000) // tiny overlap
        assert !scnvs.annotateCaller('MrBoo', cnv,  'ed')
        assert cnv.ed.best == null
        assert cnv.ed.supporting.size() == 0
        assert cnv.ed.all.size() == 0
    }  
}
