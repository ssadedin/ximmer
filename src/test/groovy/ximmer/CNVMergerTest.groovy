package ximmer

import static org.junit.Assert.*
import org.junit.Test

import gngs.*

class CNVMergerTest {

    @Test
    public void 'merge simple combination'() {
        List regions = [
            0..100,
            40..150
        ].collect { new Region('chr1', it) }
        
        CNVMerger cm = new CNVMerger(null, 0.5d)
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
        
        CNVMerger cm = new CNVMerger(null, 0.5d)
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
        
        CNVMerger cm = new CNVMerger(null, 0.5d)
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
        
        CNVMerger cm = new CNVMerger(null, 0.5d)
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
        
        CNVMerger cm = new CNVMerger(cnvs, 0.5d)        
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
        
        CNVMerger cm = new CNVMerger(cnvs, 0.5d)        
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
        
        CNVMerger cm = new CNVMerger(cnvs, 0.5d)        
        Region span = new Regions(cnvs).reduce()[0]
        List merged = cm.mergeCluster(span)
        
        // Desired outcome: 3 clusters
        assert merged.size() == 3
        assert merged[0].range == 0..10000
        assert merged[1].range == 30..150
        assert merged[2].range == 9900..9950
    }     
}
