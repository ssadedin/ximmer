package ximmer

import graxxia.*

import gngs.ToolBase
import groovy.util.logging.Log

@Log
class NormalizeTargetCoverage extends ToolBase {

    @Override
    public void run() {
        
        String cov_file = opts.i
        
        log.info "Loading $cov_file"
        
        def m = Matrix.load(cov_file, r: true)
        m = m.grep { !ean.contains('NTC') }
        
        int intervalsToSample = opts.n.toInteger()
        
        Matrix subsetted = m.transformRows { row -> 
            def subset = row[0..<intervalsToSample];  
            double mean = Stats.mean(subset); 
            return subset.collect { it / mean }
        }
        
        if(opts.b) {
            log.info "Saving to binary compressed form $opts.b"
            subsetted.saveBinary(opts.b)
        }
        else {
            System.out.withWriter { w ->
                log.info "Printing uncompressed form to stdout ..."
                subsetted.save(w)
            }
        }
    }
    
    static void main(String[] args) {
        cli('NormalizeTargetCoverage -i <Interval Summary> [-n <intervals>] [-b <binary>]', 'Read GATK interval summary and divide each sample by its mean output', args) {
            i 'GATK interval summary file', args:1, required: true
            n 'Number of intervals to sample', args:1, required: true
            b 'Binary file to save for accelerated loading', args:1, required: false
        }
    }
}
