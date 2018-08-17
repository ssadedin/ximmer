package ximmer

import gngs.Cli
import gngs.Utils
import groovy.json.JsonSlurper
import groovy.util.logging.Log

/**
 * A tool that analyses correlation output (see MultiCov) to identify controls where all samples
 * have correlation is smaller than a set threshold.
 * 
 * @author Simon Sadedin
 */
@Log
class FilterControls {

    Map correlations = null

    List<String> controls
    List<String> nonControls

    double threshold

    FilterControls(String correlationFile, List<String> controls, double threshold) {
        this.correlations =  new JsonSlurper().parseText(new File(correlationFile).readLines()[1..-1].join('\n'))
        this.controls = controls
        this.threshold = threshold
        this.nonControls = this.correlations*.key - controls
    }

    void run() {

        log.info "Total sample set is ${correlations*.key}"
        log.info "Filtering ${controls.size()} controls to correlation ${threshold}"

        log.info "${nonControls.size()} non-control samples are: $nonControls"
        
        List<String> filteredControls = controls.grep {
            isAcceptableControl(it)
        } 
        
        log.info "Selected controls are: ${filteredControls}"
        
        println nonControls.join('\n') + '\n' + filteredControls.join('\n') + '\n'
    }
    
    /**
     * Test if a control has better than the given threshold correlation to at least one test sample
     * 
     * @param control
     * @return
     */
    boolean isAcceptableControl(String control) {
        int numTestSamples = nonControls.count {  String testSample ->
            this.correlations[testSample][control] > this.threshold
        }
        
        if(numTestSamples>0) {
            log.info "Control $control is well matched to ${numTestSamples} test samples"
            return true
        }
        else {
            log.info "Control $control is not well matched to any test samples: this control will be rejected"
            return false
        }
    }

    static void main(String [] args) {

        Utils.configureSimpleLogging()

        Cli cli = new Cli(usage: 'FilterControls <options>')
        cli.with {
            corr 'JSON file containing samplewise correlations', args:1, required: true
            control 'Specify a control sample', args: Cli.UNLIMITED, required: true
            thresh 'Minimum correlation threshold to accept a control (0.9)', args: 1, required: false
        }

        OptionAccessor opts = cli.parse(args)
        if(!opts)
            System.exit(1)

        double threshold = opts.thresh? opts.thresh.toDouble() : 0.9
        new FilterControls(opts.corr, opts.controls, threshold).run()
    }
}
