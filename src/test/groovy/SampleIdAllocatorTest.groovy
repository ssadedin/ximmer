import static org.junit.Assert.*

import org.junit.Test

class SampleIdAllocatorTest {

    @Test
    public void test() {
        SampleIdAllocator sia = SampleIdAllocator.instance
        
        assert sia.newSampleId('S1') == 'XS001'
        assert sia.newSampleId('S2') == 'XS002'
        assert sia.newSampleId('S3') == 'XS003'
        
    }

}
