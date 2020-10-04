package ximmer

import gngs.Region
import groovy.transform.CompileStatic

@CompileStatic
interface OverlapCriteria {
    boolean overlaps(List<Region> r1, List<Region> r2)
}

