package ximmer

import gngs.Region
import gngs.Regions
import groovy.transform.CompileStatic

/**
 * Overlap criteria that considers two CNV calls to be the same if they
 * share at least a given fraction of their target regions, irrespective
 * of how large or small the regions are.
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class OverlapByFractionOfTargetRegions implements OverlapCriteria {
    
    Regions targetRegions
    
    double minimumFraction = 0.5

    @Override
    public boolean overlaps(Region r1, Region r2) {
        
        List<Region> r1Overlaps = targetRegions.getOverlapRegions(r1)
        List<Region> r2Overlaps = targetRegions.getOverlapRegions(r2)
        
        int mutual = (int)r1Overlaps.count { Region r -> r2Overlaps.any { r.overlaps(it) } }
        int total = (
            r1Overlaps.size() + r2Overlaps.size()
            - mutual /* the mutual ones will have been counted twice otherwise */)
        
        if(total == 0)
            return false
            
        return (mutual / total) >= minimumFraction
    }
}


