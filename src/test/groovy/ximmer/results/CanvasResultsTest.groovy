package ximmer.results

import static org.junit.Assert.*

import org.junit.Test

import gngs.*

class CanvasResultsTest {

   @Test
    public void test() {
        CanvasResults cnvs = new CanvasResults('src/test/data/test.canvas.vcf')
        Region del = new Region("chr1:30028119-72764887")
        
        assert cnvs.getOverlaps(del).size() == 1
        
        Region canvasDup = cnvs.find { it.overlaps(del) }
        assert canvasDup.type == 'DEL' // Note in CANVAS we are treating LOH as DEL
        assert canvasDup.quality > 30 && canvasDup.quality < 50
        
    }
}
