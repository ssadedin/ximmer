package ximmer

import java.util.List

import gngs.Region
import gngs.Regions
import groovy.transform.CompileStatic

class Exclusions {
    
    final static long MAX_REGION_WAIT_TIME_MS = 300000
    
    Regions targetRegions
    
    /**
     * The current set of regions that is excluded
     */
    Regions regions
    
    List<Region> pendingReservations = []
   
    public Exclusions(Regions targetRegions, Regions regions) {
        super();
        this.targetRegions = targetRegions;
        this.regions = regions;
    }

    /**
     * Check if the given region can be accepted as a simulated CNV
     * region, by first checking if it overlaps any existing regions
     * and then also passing it to the given evaluator function which 
     * can expand the region.
     * 
     * @param region    Region to check
     * @param evaluator Evaluation function to further process the region
     * @return  null if not accepted, or the post-evaluation region if accepted
     */
    @CompileStatic
    Region tryReserve(Region region, Closure evaluator) {
        
        long startMs = System.currentTimeMillis()
        
        Region paddedRegion = padRegion(region)
        
        synchronized(regions) {
            
            if(regions.overlaps(paddedRegion))
                return null
                
            List overlapping = pendingReservations.grep { Region r -> r.overlaps(paddedRegion) }
            if(!overlapping.isEmpty())
                return null
                
            pendingReservations << paddedRegion
        }
        
        Region acceptedRegion = (Region)evaluator(region)
        synchronized(regions) {
            pendingReservations.remove(paddedRegion)
            if(acceptedRegion != null)
                regions.addRegion(padRegion(acceptedRegion))
        }
        return acceptedRegion
    }

    /**
     * Add a target region upstream  and downstream of the given region and
     * return the result
     */
    private Region padRegion(Region region) {
        IntRange previous = targetRegions.previousRange(region.chr, region.from)
        IntRange next = targetRegions.nextRange(region.chr, region.to)

        Region paddedRegion = new Region(
                region.chr,
                previous ? previous.from : region.from,
                next ? next.to : region.to)
        return paddedRegion
    }
    
}
