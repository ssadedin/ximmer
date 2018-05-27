package ximmer

import static org.junit.Assert.*

import gngs.Region
import gngs.Regions
import org.junit.Test

class ExclusionsTest {

    Regions targetRegions = [
        r('chr1:100-200'),
        r('chr1:300-400'),
        r('chr1:600-700'),
        r('chr1:900-1000'),
    ] as Regions
        
    @Test
    public void testNoExclusions() {
        
        Exclusions exc = new Exclusions(new Regions(), new Regions())
        
        def a = exc.tryReserve(r('chr1:100-200'), { it }) 
        assert a != null
        assert a.chr == 'chr1'
        assert a.range == 100..200
    }
    
    @Test
    public void testSimpleExclusion() {
        Exclusions exc = new Exclusions(new Regions(), new Regions([r('chr1:100-200')]))
        def a = exc.tryReserve(r('chr1:100-200'), { it }) 
        assert a == null
    }
   
    @Test
    public void 'region should be expanded 1 target upstream and downstream'() {
        Exclusions exc = new Exclusions(targetRegions, new Regions())
        def a = exc.tryReserve(r('chr1:300-400'), { it })
        assert a != null
        assert exc.regions[0].from == 100
        assert exc.regions[0].to == 700
    }
    
    @Test
    public void 'evaluator expanded region'() {
        Exclusions exc = new Exclusions(new Regions(), new Regions())
        
        def a = exc.tryReserve(r('chr1:300-400'), { 
            new Region('chr1:100-700')
        })         
        
        assert a != null
        assert exc.regions[0].from == 100
        assert exc.regions[0].to == 700        
    }
    
    @Test
    public void 'padding is added to evaluator expansion'() {
        Exclusions exc = new Exclusions(targetRegions, new Regions())
        
        def a = exc.tryReserve(r('chr1:300-400'), { 
            new Region('chr1:100-700')
        })         
        
        assert a != null
        assert exc.regions[0].from == 100
        assert exc.regions[0].to == 1000
    } 
    
    Region r(String value) {
        return new Region(value)
    }

}
