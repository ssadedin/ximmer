package ximmer.results

import static org.junit.Assert.*

import org.junit.Test

import gngs.*

class CanvasResultsTest {

   @Test
    public void test() {
        CanvasResults cnvs = new CanvasResults('src/test/data/test.canvas.vcf')
//        Region del = new Region("chr1:30028119-72764887")
        Region del = new Region("chr1:144025582-144031152")
        
        assert cnvs.getOverlaps(del).size() == 1
        
        Region canvasDup = cnvs.find { it.overlaps(del) }
        assert canvasDup.type == 'DEL' // Note in CANVAS we are treating LOH as DEL
        assert canvasDup.quality > 1.0 && canvasDup.quality < 3.0
        assert canvasDup.cn == 1
        
    }
}
