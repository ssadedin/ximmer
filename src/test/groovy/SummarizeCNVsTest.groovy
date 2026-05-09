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
        
        // Set up fake caller data the way annotateCaller would
        cnv.ed = [best: null, supporting: [], all: []]
        cnv.xhmm = [best: null, supporting: [], all: []]
        
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
        
        // cnvToMap expects cnv[caller] to be a Map, as set by annotateCaller
        cnv.ed = [best: ed, supporting: [ed], all: [ed]]
        cnv.xhmm = [best: null, supporting: [], all: []]
        
        Map data = scnvs.cnvToMap(callers, [], scnvs.computeColumns(callers,[]), cnv)
        
        println data
        
        Utils.table([data])
        
        assert data.ed == 'TRUE'
        assert data.xhmm == 'FALSE'
        assert data.calls.ed[0] == [900,1300,10]
    }
    
    @Test
    void 'test single caller are annotated correctly'() {
        Region edCall = new Region('chr1', 1000..2000, quality: 100, sample: 'MrBoo')
        edCall.caller = 'ed'
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // annotateCaller now takes (Region, String) — sample from cnv.sample
        Region cnv = new Region('chr1', 1200..1800)
        cnv.sample = 'MrBoo'
        cnv.cnvs = [edCall] as Set
        
        assert scnvs.annotateCaller(cnv, 'ed')
        assert cnv.ed.best.is(edCall)
        assert cnv.ed.supporting.size() == 1
        assert cnv.ed.all.size() == 1
    }
    
    @Test
    void 'test single call only supported when mutal overlap'() {
        Region edCall = new Region('chr1', 1200..1250, quality: 100, sample: 'MrBoo')
        edCall.caller = 'ed'
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // annotateCaller now takes (Region, String)
        Region cnv = new Region('chr1', 1000..2000) // tiny overlap
        cnv.sample = 'MrBoo'
        cnv.cnvs = [] as Set  // no calls in the merged cluster with caller='ed'
        
        assert !scnvs.annotateCaller(cnv, 'ed')
        assert cnv.ed.best == null
        assert cnv.ed.supporting.size() == 0
        assert cnv.ed.all.size() == 1
    } 
    
    @Test
    void 'test single call only supported when correct sample'() {
        Region edCall = new Region('chr1', 1200..1250, quality: 100, sample: 'MsFoo')
        edCall.caller = 'ed'
        scnvs.results = [
            'ed' : new Regions([edCall])
        ]
        
        // annotateCaller now takes (Region, String)
        Region cnv = new Region('chr1', 1000..2000) // tiny overlap
        cnv.sample = 'MrBoo'  // different sample!
        cnv.cnvs = [] as Set
        
        assert !scnvs.annotateCaller(cnv, 'ed')
        assert cnv.ed.best == null
        assert cnv.ed.supporting.size() == 0
        assert cnv.ed.all.size() == 0  // no calls from this caller for this sample
    }  
}
