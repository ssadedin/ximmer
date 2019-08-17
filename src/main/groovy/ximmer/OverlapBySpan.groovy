package ximmer

import gngs.Region
import groovy.transform.CompileStatic

@CompileStatic
class OverlapBySpan implements OverlapCriteria {

    double minimumFraction = 0.5
    
    @Override
    public boolean overlaps(Region r1, Region r2) {
        return r1.mutualOverlap(r2) > minimumFraction;
    }
}

