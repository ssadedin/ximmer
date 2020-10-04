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
    Regions mergeCluster(final Region cluster) {
        
        List<Region> overlaps = cnvs.getOverlapRegions(cluster)
        
        // For each combination of overlapping region, determine mutual overlap
        Regions finalMerge = overlaps as Regions
        while(true) {
            Regions newMerge = mergeMutualOverlapping(finalMerge)
            if(newMerge.numberOfRanges == finalMerge.numberOfRanges)
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
    @CompileStatic
    Regions mergeMutualOverlapping(final Regions overlaps) {
        Regions combined = new Regions()
        final Set allMerged = new HashSet<Region>()
        for(Region cnv1 : overlaps) {
            Region newCluster = findCNVCluster(cnv1, overlaps)
            
            boolean wasMerged = false
            combined = combined.collect { Region existingCluster ->
                boolean overlappingCluster = 
                    this.overlapCriteria.overlaps([existingCluster],[newCluster])
                if(overlappingCluster) {
                    wasMerged = true
                    return mergeClusters(existingCluster, newCluster)
                }
                else
                    return existingCluster
            } as Regions

            if(!wasMerged) {
                assert newCluster != null
                combined.addRegion(newCluster)
            }
        }
        
        return combined.uniquify()
    }

    @CompileStatic
    private Region mergeClusters(Region existingCluster, Region newCluster) {
        Region merged = existingCluster.union(newCluster)
        HashSet cnvs = (HashSet)existingCluster['cnvs']
        HashSet otherCnvs = (HashSet)newCluster['cnvs']
        HashSet mergedCnvs = new HashSet()
        if(cnvs)
            mergedCnvs.addAll(cnvs)
        if(otherCnvs)
            mergedCnvs.addAll(otherCnvs)
        merged['cnvs'] = mergedCnvs
        return merged
    }

    /**
     * Calculate a region that represents the given CNV merged of all the given regions that satisfy 
     * the overlap criteria with it
     * 
     * @param cnv1
     * @param overlaps
     * @return
     */
    private Region findCNVCluster(Region cnv1, Regions overlaps) {
        
        String cnvCaller = cnv1.properties.getOrDefault('caller', 'all').toString()
        Map<String,List<Region>> otherCallerCalls = overlaps.groupBy { 
            it.properties.getOrDefault('caller','all').toString()
        }
        
        List<Region> cnvCallerCalls = otherCallerCalls[cnvCaller]
        
        Region newCluster = cnv1
        for(Map.Entry otherCallerEntry : otherCallerCalls) {
            String otherCaller = otherCallerEntry.key
            List<Region> callerCNVs = otherCallerEntry.value

            // We are only interested in merging across different groups, not the same
            if(otherCaller == cnvCaller)
                continue
            
            if(this.overlapCriteria.overlaps(cnvCallerCalls, callerCNVs)) {

                Regions newClusterRegions = 
                    new Regions(callerCNVs)
                        .addRegion(cnv1)
                        .reduce()

                newCluster = newClusterRegions.getSpan(cnv1.chr)
                HashSet allCnvs = new HashSet()
                for(Region cnv2 : cnvCallerCalls) {
                    allCnvs.addAll(cnv2)
                }

                for(Region cnv2 : callerCNVs) {
                    allCnvs.addAll(cnv2.cnvs?:[cnv2])
                }
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
