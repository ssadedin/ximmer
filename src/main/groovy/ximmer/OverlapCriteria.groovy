package ximmer

import gngs.Region
import groovy.transform.CompileStatic

/**
 * Evaluates the overlap between two regions via different methods
 */
@CompileStatic
interface OverlapCriteria {

    boolean overlaps(Region r1, Region r2)
    
    double calculateOverlap(Region r1, Region r2)

}

