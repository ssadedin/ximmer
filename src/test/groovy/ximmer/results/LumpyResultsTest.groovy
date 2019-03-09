package ximmer.results

import static org.junit.Assert.*

import org.junit.Test

import gngs.*

class LumpyResultsTest {


    @Test
    public void test() {
        String vcf_path = 'src/test/data/test.small.lumpy.vcf'
        LumpyResults lr = new LumpyResults('src/test/data/test.small.lumpy.vcf')
        Region dup = new Region("chr2:7283-7743")

        assert lr.getOverlaps(dup).size() == 1

        Region lumpyDup = lr.find { it.overlaps(dup) }
        assert lumpyDup.type == 'DUP'
        assert lumpyDup.quality > 9
        assert lumpyDup.size() == 368

    }

    @Test
    void 'correct sample id inferred'() {
        LumpyResults lr = new LumpyResults('src/test/data/test.small.lumpy.vcf')
        assert lr[0].sample == 'TESTSAMPLE' // See VCF
    }

}
