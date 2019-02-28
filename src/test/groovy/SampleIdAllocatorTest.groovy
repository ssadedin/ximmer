import static org.junit.Assert.*

import org.junit.Test

class SampleIdAllocatorTest {

    @Test
    public void test() {
        SampleIdAllocator sia = SampleIdAllocator.instance
        
        assert sia.newSampleId('S1') == 'S1'
        assert sia.newSampleId('S2') == 'S2'
        assert sia.newSampleId('S3') == 'S3'
        
    }

}
