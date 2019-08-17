package ximmer

import gngs.*
import groovy.transform.CompileStatic


@CompileStatic
interface OverlapCriteria {
    boolean overlaps(Region r1, Region r2)
}

@CompileStatic
class OverlapBySpan implements OverlapCriteria {

    double minimumFraction = 0.5
    
    @Override
    public boolean overlaps(Region r1, Region r2) {
        return r1.mutualOverlap(r2) > minimumFraction;
    }
}


@CompileStatic
class OverlapByFractionOfTargetRegions implements OverlapCriteria {
    
    Regions targetRegions
    
    double minimumFraction = 0.5

    @Override
    public boolean overlaps(Region r1, Region r2) {
        
        List<Region> r1Overlaps = targetRegions.getOverlapRegions(r1)
        List<Region> r2Overlaps = targetRegions.getOverlapRegions(r2)
        
        int mutual = (int)r1Overlaps.count { Region r -> r in r2Overlaps }
        
        int total = (
            r1Overlaps.size() + r2Overlaps.size()
            - mutual /* the mutual ones will have been counted twice otherwise */)
        
        return (mutual / total) >= minimumFraction
    }
}

class CNVMerger {
    
    final Regions cnvs
    
    OverlapCriteria overlapCriteria
    
    public CNVMerger(Regions cnvs, OverlapCriteria overlapCriteria) {
        super();
        this.cnvs = cnvs;
        this.overlapCriteria = overlapCriteria
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
     * Combine all pairs regions in the list that satisfy the {@link #overlapCriteria}
     * to overlap with each other.
     * <p>
     * Note that the same region may be included multiple times in the result.
     * 
     * @param overlapping   List of regions to combine
     * @return  List of combined regions, or empty list if none could be combined
     */
    List<Region> mergeMutualOverlapping(final List<Region> overlaps) {
        List combined = []
        final Set allMerged = new HashSet<Region>()
        for(int i=0; i<overlaps.size(); ++i) {
            Region cnv1 = overlaps[i]
            Region newCluster = findCNVCluster(cnv1, overlaps)
            
            boolean wasMerged = false
            combined = combined.collect { Region existingCluster ->
                if(this.overlapCriteria.overlaps(existingCluster,newCluster)) {
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

    /**
     * Calculate a region that represents the given CNV merged of all the given regions that satisfy 
     * the overlap criteria with it
     * 
     * @param cnv1
     * @param overlaps
     * @return
     */
    private Region findCNVCluster(Region cnv1, List overlaps) {
        Region newCluster = cnv1
        for(int j=0; j< overlaps.size(); ++j) {
            Region cnv2 = overlaps[j]
            if(cnv1.is(cnv2))
                continue
            
            if(this.overlapCriteria.overlaps(cnv1, cnv2)) {
                newCluster = newCluster.union(cnv2)
                HashSet allCnvs = new HashSet()
                allCnvs.addAll(cnv1.cnvs?:[cnv1])
                allCnvs.addAll(cnv2.cnvs?:[cnv2])
                newCluster.cnvs = allCnvs
            }
        }
        return newCluster
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
