package ximmer

import gngs.Region
import groovy.transform.CompileStatic

@CompileStatic
class OverlapBySpan implements OverlapCriteria {

    double minimumFraction = 0.5
    
    @Override
    public boolean overlaps(Region r1, Region r2) {
        final int maxBreakEndDistance = 3i;
        
        if((r1['type'] == 'BND') ^ (r2['type'] == 'BND')) {
            return Math.abs(r1.from - r2.from) < maxBreakEndDistance || Math.abs(r2.to - r1.to) < maxBreakEndDistance
        }

        return r1.mutualOverlap(r2) > minimumFraction;
    }

    @Override
    public double calculateOverlap(Region r1, Region r2) {
        
        if((r1['type'] == 'BND') ^ (r2['type'] == 'BND')) {
            return 0d
        }
        
        return r1.mutualOverlap(r2) > minimumFraction;
    }
}

