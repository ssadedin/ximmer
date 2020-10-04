package ximmer

import gngs.Region
import gngs.Regions
import groovy.transform.CompileStatic

@CompileStatic
class OverlapBySpan implements OverlapCriteria {

    double minimumFraction = 0.5
    
    @Override
    public boolean overlaps(List<Region> r1, final List<Region> r2) {
        Regions r1r = r1 as Regions
        Regions r2r = r2 as Regions
        
        Regions ix = r1r.intersect(r2r)
        
        double totalFracOverlap = Math.min((double)(ix.size() / r1r.size()), (double)(ix.size() / r2r.size()))

        return totalFracOverlap > minimumFraction;
    }
}

