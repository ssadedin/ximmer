package ximmer

import gngs.Region
import groovy.transform.CompileStatic

@CompileStatic
interface OverlapCriteria {
    boolean overlaps(Region r1, Region r2)
}

