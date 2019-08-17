package ximmer

import static org.junit.Assert.*
import org.junit.Test

import gngs.*

class CNVMergerTest {
    
    OverlapBySpan spanOverlap = new OverlapBySpan(minimumFraction: 0.5)
    
    @Test
    public void 'merge simple combination'() {
        List regions = [
            0..100,
            40..150
        ].collect { new Region('chr1', it) }
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        List merged = cm.mergeMutualOverlapping(regions)
        assert merged.size() == 1
        assert merged[0].from == 0
        assert merged[0].to == 150
    }
    
    @Test
    public void 'merge two with one separate'() {
        List regions = [
            0..100,
            40..150,
            145..150
        ].collect { new Region('chr1', it) }
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        List merged = cm.mergeMutualOverlapping(regions)
        assert merged.size() == 2
        assert merged[0].from == 0
        assert merged[0].to == 150
        assert merged[1].range == 145..150
    }
    
    @Test
    public void 'large and tiny do not merge'() {
        List regions = [
            0..100,
            50..54,
        ].collect { new Region('chr1', it) }
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        List merged = cm.mergeMutualOverlapping(regions)
        assert merged.size() == 2
        assert merged[0].range == 0..100
        assert merged[1].range == 50..54
    } 
    
    @Test
    public void 'merge of cascading regions'() {
        List regions = [
            0..100,
            30..130,
            70..170,
            110..210
        ].collect { new Region('chr1', it) }
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        List merged = cm.mergeMutualOverlapping(regions)
        assert merged.size() == 1
        assert merged[0].range == 0..210
    } 
     
    @Test
    void 'simple cluster of two merges to one'() {
        Regions cnvs = [
            0..100,
            30..130,
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        List merged = cm.mergeCluster(new Region('chr1:0-130'))
        
        assert merged.size() == 1
        assert merged[0].range == 0..130
    }
    
    @Test
    void 'cluster of three merges to one'() {
        Regions cnvs = [
            0..100,
            30..130,
            50..150,
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        List merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.size() == 1
        assert merged[0].range == 0..150
    }
    
    @Test
    void 'complex merge of many regions'() {
        Regions cnvs = [
            // large overlapping call
            0..10000,
            
            // one small cluster at the start
            30..130, 
            50..150,
            60..115,
            
            // second small cluster at the end
            9900..9940,
            9920..9950,
            9930..9950
            
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        Region span = new Regions(cnvs).reduce()[0]
        List merged = cm.mergeCluster(span)
        
        // Desired outcome: 3 clusters
        assert merged.size() == 3
        assert merged[0].range == 0..10000
        assert merged[1].range == 30..150
        assert merged[2].range == 9900..9950
    }     
    
    OverlapByFractionOfTargetRegions targetOverlap = new OverlapByFractionOfTargetRegions(
        minimumFraction: 0.5, targetRegions:
        [
        20..40,
        60..80,
        100..120,
        140..160
    ].collect { new Region("chr1", it) } as Regions)

    @Test
    void 'cluster of two merges to one by target'() {
        Regions cnvs = [
            0..120, // overlaps first 3 targets
            60..160, // overlaps last 3 targets
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        List merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.size() == 1
        assert merged[0].range == 0..160
    }
    
    @Test
    void 'one third of target regions does not merge by target'() {
        Regions cnvs = [
            0..40, // overlaps first target
            0..120, // overlaps last 3 targets
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        List merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.size() == 2
    }
    
    @Test
    void 'one large and two disjoint overlap by target'() {
        Regions cnvs = [
            0..80, // overlaps first 2 targets
            100..160, // overlaps last 2 targets
            0..160, // overlaps all 4 targets
        ].collect { new Region('chr1', it) } as Regions
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        List merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.size() == 1
    }
}
