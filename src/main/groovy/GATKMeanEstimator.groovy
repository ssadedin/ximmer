import gngs.GATKIntervalSummary
import groovy.util.logging.Log;

/**
 * Estimates sample means by reading interval summary file output by GATK
 * 
 * @author simon
 */
@Log
class GATKMeanEstimator {
    
    File dir = null
    
    Map<String,File> intervalFiles = [:]
    
    GATKMeanEstimator(String dir) {
        this.dir = new File(dir)
        
        this.intervalFiles = this.dir.listFiles().grep { 
            it.name.endsWith(".sample_interval_summary") 
        }.collectEntries {  File f ->
            // Extract the sample name from the first line
            String firstLine = f.withReader { r -> r.readLine() }
            
            // Find the column that ends with "_granular_median"
            String meanColumn = firstLine.tokenize("\t").find { it.endsWith("_granular_median") }
            
            String sample = meanColumn.replaceAll('_granular_median','')
            
            [ sample, f]
        }
    }
    
    Map<String,Double> calculateMeans(List<String> samples) {
        samples.collectEntries { sample ->
            GATKIntervalSummary summary = new GATKIntervalSummary(intervalFiles[sample].absolutePath)
            summary.load()
            
            log.info "Loading mean for $sample"
            [
                sample,
                summary*.average_coverage.sum() / summary.numberOfRanges
            ]
        }
    }
}
