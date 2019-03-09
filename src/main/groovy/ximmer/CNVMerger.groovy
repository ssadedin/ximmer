package ximmer

import gngs.*
import groovy.transform.CompileStatic

class CNVMerger {
    
    final Regions cnvs
    
    final double fracOverlap
    
    public CNVMerger(Regions cnvs, double fracOverlap) {
        super();
        this.cnvs = cnvs;
        this.fracOverlap = fracOverlap;
    }

    @CompileStatic
    Regions merge() {
        List<Region> result = new ArrayList(cnvs.numberOfRanges)
        Regions reduced = cnvs.reduce()
        for(Region cluster in reduced) {
            result.addAll(mergeCluster(cluster))
        }
        return result as Regions
    }
    
    @CompileStatic
    List<Region> mergeCluster(Region cluster) {
        
        List<Region> overlaps = cnvs.getOverlapRegions(cluster)
        
        // For each combination of overlapping region, determine mutual overlap
        List<Region> finalMerge = overlaps
        while(true) {
            List<Region> newMerge = mergeMutualOverlapping(finalMerge)
            if(newMerge.size() == finalMerge.size())
                break
            finalMerge = newMerge
        }
        return finalMerge 
    }
        
    /**
     * Combine all pairs regions in the list that have at least {@link #fracOverlap} mutual
     * overlap with each other.
     * <p>
     * Note that the same region may be included multiple times in the result.
     * 
     * @param overlapping   List of regions to combine
     * @return  List of combined regions, or null if none could be combined
     */
    List<Region> mergeMutualOverlapping(final List<Region> overlaps) {
        List combined = []
        final Set allMerged = new HashSet<Region>()
        for(int i=0; i<overlaps.size(); ++i) {
            Region cnv1 = overlaps[i]
            Region newCluster = cnv1
            for(int j=0; j< overlaps.size(); ++j) {
                if(i==j)
                    continue
                Region cnv2 = overlaps[j]
                if(cnv1.mutualOverlap(cnv2) > fracOverlap) {
                    newCluster = newCluster.union(cnv2)
                    HashSet allCnvs = new HashSet()
                    allCnvs.addAll(cnv1.cnvs?:[cnv1])
                    allCnvs.addAll(cnv2.cnvs?:[cnv2])
                    newCluster.cnvs = allCnvs 
                }
            }
            
            boolean wasMerged = false
            combined = combined.collect { Region existingCluster ->
                if(existingCluster.mutualOverlap(newCluster) > fracOverlap) {
                    wasMerged = true
                    return existingCluster.union(newCluster)
                }
                else
                    return existingCluster
            }
            if(!wasMerged) {
                assert newCluster != null
                combined.add(newCluster)
            }
        }
        
        return combined
    }
        
        
//        // Algorithm: iterate through the cnvs, accumulating the set of merged regions as we go,
//        // and marking each CNV that we merge with a 'merged' flag
//        //
//        // If we encounter a 'merged' cnv later as we iterate, we skip it because it has already
//        // been merged.
//        for(Region cnv in cnvs) {
//            List<Region> overlaps = cnvs.getOverlapRegions().grep { Region r -> r.mutualOverlap(cnv) > fracOverlap }
//            
//            Region existingMerge = overlaps.find { it.merged }
//            Region cnvMerge
//            if(existingMerge) {
//                cnvMerge = existingMerge.union(cnv)
//                cnv.merged = cnvMerge
//            }
//            else {
//                cnvMerge = cnv
//                for(Region overlap in overlaps) {
//                    cnvMerge = cnvMerge.union(overlaps)
//                }
//            }
//        }        
}
