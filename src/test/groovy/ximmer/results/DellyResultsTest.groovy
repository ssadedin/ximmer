package ximmer.results

import static org.junit.Assert.*

import gngs.Region
import org.junit.Test

class DellyResultsTest {

    @Test
    public void test() {
        DellyResults dr = new DellyResults('src/test/data/test.delly.vcf')
        Region dup = new Region("chr1:10200-10443")
        
        assert dr.getOverlaps(dup).size() == 1
        
        Region dellyDup = dr.find { it.overlaps(dup) }
        assert dellyDup.type == 'DUP'
        assert dellyDup.quality > 10
        
    }
    
    @Test
    void 'test no inversions'() {
        DellyResults dr = new DellyResults('src/test/data/test.delly.vcf')
        
        assert !dr*.type.any { it.contains('INV') }
        
        Region inv = new Region("chr1:782758-782858")
        
        assert dr.getOverlaps(inv).isEmpty()
    }
    
    @Test
    void 'correct sample id inferred'() {
        DellyResults dr = new DellyResults('src/test/data/test.delly.vcf')
        assert dr[0].sample == 'TESTSAMPLE' // See VCF
    }
}
