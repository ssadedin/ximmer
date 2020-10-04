package ximmer

import static org.junit.Assert.*
import org.junit.Test

import gngs.*

class CNVMergerTest {
    
    OverlapBySpan spanOverlap = new OverlapBySpan(minimumFraction: 0.5)
    
    @Test
    public void 'merge simple combination'() {
        Regions regions = cnvRegions(
            0..100,
            40..150
        )
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        Regions merged = cm.mergeMutualOverlapping(regions)
        assert merged.numberOfRanges == 1
        assert merged[0].from == 0
        assert merged[0].to == 150
    }
    
    @Test
    public void 'merge two with one separate'() {
        Regions regions = cnvRegions(
            0..100,
            40..150,
            145..150
        )
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        Regions merged = cm.mergeMutualOverlapping(regions)
        assert merged.numberOfRanges == 2
        assert merged[0].from == 0
        assert merged[0].to == 150
        assert merged[1].range == 145..150
    }
    
    @Test
    public void 'large and tiny do not merge'() {
        Regions regions = cnvRegions(
            0..100,
            50..54,
        )
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        Regions merged = cm.mergeMutualOverlapping(regions)
        assert merged.numberOfRanges == 2
        assert merged[0].range == 0..100
        assert merged[1].range == 50..54
    } 
    
    @Test
    public void 'merge of cascading regions'() {
        Regions regions = cnvRegions(
            0..100,
            30..130,
            70..170,
            110..210
        )
        
        CNVMerger cm = new CNVMerger(null, spanOverlap)
        Regions merged = cm.mergeMutualOverlapping(regions)
        assert merged.numberOfRanges == 1
        assert merged[0].range == 0..210
    } 
     
    @Test
    void 'simple cluster of two merges to one'() {
        Regions cnvs = cnvRegions(
            0..100,
            30..130,
        )
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:0-130'))
        
        assert merged.numberOfRanges == 1
        assert merged[0].range == 0..130
    }
    
    @Test
    void 'cluster of three merges to one'() {
        Regions cnvs = cnvRegions(
            0..100,
            30..130,
            50..150,
        )
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.numberOfRanges == 1
        assert merged[0].range == 0..150
    }
    
    @Test
    void 'complex merge of many regions'() {
        Regions cnvs = cnvRegions(
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
        )
        
        CNVMerger cm = new CNVMerger(cnvs, spanOverlap)        
        Region span = new Regions(cnvs).reduce()[0]
        Regions merged = cm.mergeCluster(span)
        
        // Desired outcome: 3 clusters
        assert merged.numberOfRanges == 3
        assert merged[0].range == 0..10000
        assert merged[1].range == 30..150
        assert merged[2].range == 9900..9950
    }     
    
    OverlapByFractionOfTargetRegions targetOverlap = new OverlapByFractionOfTargetRegions(
        minimumFraction: 0.5, targetRegions:
        [ 20..40, 60..80, 100..120, 140..160, 180..200, ]
        .collect { new Region("chr1", it) } as Regions)

    @Test
    void 'cluster of two merges to one by target'() {
        Regions cnvs = cnvRegions(
            0..120, // overlaps first 3 targets
            60..160, // overlaps last 3 targets
        )
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.numberOfRanges == 1
        assert merged[0].range == 0..160
    }
    
    @Test
    void 'one third of target regions does not merge by target'() {
        Regions cnvs = cnvRegions(
            0..40, // overlaps first target
            0..120, // overlaps last 3 targets
        )
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.numberOfRanges == 2
    }
    
    @Test
    void 'one large and two disjoint overlap by target'() {
        Regions cnvs = cnvRegions(
            0..80, // overlaps first 2 targets
            100..160, // overlaps last 2 targets
            0..160, // overlaps all 4 targets
        )
        
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:0-150'))
        
        assert merged.numberOfRanges == 1
    }

    /**
     * When several small calls by one caller overlap a large
     * one by a different caller, the calls should be combined
     * to test the overlap criteria:
     * 
     *    |-ed-| |-ed-| |-ed--|
     *    |--------xhmm--------|
     *             ⬇
     *    |--------------------|
     */
    @Test
    void 'two disjoint by one caller and one large'() {
        // [ 20..40, 60..80, 100..120, 140..160, 180..200, ]

        Regions edCnvs = cnvRegions(
            20..80, // overlaps first 2 targets
            140..200, // overlaps last 2 targets
            caller:'ed'
        )
       
        Regions xhmmCnvs = cnvRegions(
            20..200, // overlaps all 4 targets
            caller: 'xhmm'
        )
        
        Regions cnvs = edCnvs + xhmmCnvs
        CNVMerger cm = new CNVMerger(cnvs, targetOverlap)        
        Regions merged = cm.mergeCluster(new Region('chr1:20-200'))
        
        assert merged.numberOfRanges == 1
        assert merged[0].from == 20
        assert merged[0].to == 200
        assert merged[0].cnvs.size() == 3
    }
    
    @Test
    void 'cnv cluster with 3 nonoverlapping stays 3'() {
        Regions clusters = cnvRegions(0..10000, 30..150, 9900..9950)
        CNVMerger cm = new CNVMerger(clusters, targetOverlap)        
        Region cnv1 = clusters[-1]
        
        def result = cm.findCNVCluster(cnv1, clusters)
        assert result.from == 9900
        assert result.to == 9950
    }
    
    Regions cnvRegions(Map attributes=[:], IntRange... cnvs) {
        String chr = attributes.getOrDefault('chr','chr1')
        int i = 1
        return cnvs.collect { new Region(chr, it) }
             .each {
                attributes.each { k,v ->
                    it[k] = v
                }
                if(!attributes['caller']) {
                    it['caller'] = "caller${i++}".toString()
                }
            } 
            .asType(Regions)
     }
}
