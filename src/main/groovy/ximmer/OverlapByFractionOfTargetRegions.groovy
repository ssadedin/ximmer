package ximmer

import gngs.Region
import gngs.Regions
import groovy.transform.CompileStatic

@CompileStatic
class OverlapByFractionOfTargetRegions implements OverlapCriteria {
    
    Regions targetRegions
    
    double minimumFraction = 0.5

    @Override
    public boolean overlaps(List<Region> r1, List<Region> r2) {
        
        List<Region> r1Overlaps = r1.collectMany { Region r ->
            def o = targetRegions.getOverlapRegions(r)
            return o
        }

        List<Region> r2Overlaps = r2.collectMany { Region r ->
            targetRegions.getOverlapRegions(r)
        }
        
        int mutual = (int)r1Overlaps.count { Region r -> r2Overlaps.any { it.from == r.from } }
        
        int total = (
            r1Overlaps.size() + r2Overlaps.size()
            - mutual /* the mutual ones will have been counted twice otherwise */)
        
        return (mutual / total) >= minimumFraction
    }
}


