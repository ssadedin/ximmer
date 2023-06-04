package ximmer

import gngs.*
import groovy.transform.CompileStatic

/**
 * Implements logic for merging multiple CNV calls (typically, from different CNV callers) into 
 * one overarching harmonised CNV call
 * <p>
 * There are a range of factors that are important to consider in this process. The fine grained
 * logic is abstracted to an {@link OverlapCriteria} interface to allow customisation.
 * <p>
 * The high level logic consists of the following process:
 * 
 * <li> Flatten the ranges of all calls to a single flattened set of non-overlapping ranges, eg:
 * <pre>
 *        |------------|
 *               |--------|
 *               ⬇
 *        |---------------|
 * </pre>
 * <li>Each flattened range (or "cluster") represents a candidate combined call. However, some such combined
 *     calls are in apppropriate, so each such call is then further processed, by the
 *     {@link #mergeCluster} method to split it into separate regions that satisfy
 *     the configured {@link #overlapCriteria}. This process is iterative: 
 *     <ol>
 *        <li>all combinatorial pairs of regions within the cluster are tested for mutual overlap
 *        <li>all mutual overlaps are merged to a single call
 *        <li>the remaining merged calls are then fed back into the step 1
 *        <li>when an iteration occurs where no merges occur, the process ends
 *     </ol>
 * For example, consider 4 cnv callers like so:
 * <pre>
 *        |--a---| |-b--|
 *        |-c-| |----d---|
 *               ⬇
 *        |-a,c--| 
 *              |--b,d---|
 * </pre>
 * In the above example, (a) and (d) are clearly not sufficiently overlapping to warrant combination. The 
 * overlap criteria would likely reject this merge, which severs the link between them and
 * causes the cluster to split into two. This turns out to be important to stop the merging from creating
 * giant artefactual clusters through chaining of coincidental overlaps.
 * 
 * @author Simon Sadedin
 */
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
        for(Region ol in overlaps) {
            ol['cnvs'] = [ol]
        }
        
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
                    Region mergedCluster = existingCluster.union(newCluster)
                    mergedCluster.cnvs = (existingCluster.cnvs + newCluster.cnvs).unique()
                    return mergedCluster
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
                newCluster = mergeCNVClusters(newCluster,cnv1,cnv2)
            }
        }
        return newCluster
    }
    
    private Region mergeCNVClusters(Region cluster, Region cnv1, Region cnv2) {
        Region newCluster
        
        // Break ends can introduce a scenario where the two CNVs are considered
        // part of the same cluster but do not actually overlap
        if(cnv1.overlaps(cnv2)) {
            newCluster = cluster.union(cnv2)
        }
        else {
            newCluster = [cnv1,cnv2].max { it.size() }
        }

        HashSet allCnvs = new HashSet()
        allCnvs.addAll(cnv1.cnvs?:[cnv1])
        allCnvs.addAll(cnv2.cnvs?:[cnv2])
        newCluster.cnvs = allCnvs        
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
