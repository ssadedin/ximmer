package ximmer

import gngs.ToolBase
import graxxia.IntegerStats
import graxxia.Matrix
import graxxia.Stats
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool

@Log
class CalculateCombinedStatistics extends ToolBase {

    @Override
    public void run() {
        List<String> intervalStats = opts.arguments().grep { it.endsWith('.sample_interval_summary') }

        List<String> stats = opts.arguments().grep { it.endsWith('.stats.tsv') }
        
        assert stats.size() == intervalStats.size() :
             'One or more files provided as a sample interval summary did not have an associated coverage file or files were ambiguous'
        
        log.info "Loaded ${intervalStats.size()} summary files ..."

        int threads = opts.threads?:2
        log.info "Loading ${intervalStats.size()} summary files using concurrency $threads..."
        Matrix allCovs = GParsPool.withPool(threads) {
             Matrix.concat(intervalStats.collectParallel { loadCovs(it) })
        }
        
        assert allCovs.sample.every { s -> stats.count { new File(it).name.startsWith(s + '.') } == 1 } :
             'One or more files provided as a sample interval summary did not have an associated coverage file or files were ambiguous'
        
        Matrix allStats = loadCoverageStats(allCovs.sample, stats)

        log.info "Loaded ${stats.size()} coverage stats files ..."
        
        log.info "Computing row correlations ..."
        Matrix correlations = allCovs.rowCorrelations
        
        log.info "Computing coefficient of variation statistics ..."
        Matrix norm = allCovs.normaliseRows().normaliseColumns()
        List<Stats> targetRegionStats = norm.columns.collect { Stats.from(it) }
        
        IntegerStats coeffvStats = new IntegerStats(100, targetRegionStats*.standardDeviation.collect { it*100})
        
        correlations.save(opts['corrTSV'], r:true)
        log.info "Saved row correlations to $opts.corrTSV"
        
        writeCorrelationJS(opts.corrJS, correlations)
        log.info "Saved row correlations in Javascript/JSON format to $opts.corrJS"
        
        writeCovsJS(opts.covJS, allStats)
        log.info "Saved combined coverage Javascript/JSON format to $opts.covJS"
        
        writeCoeffVStats(opts.coeffvJS, coeffvStats)
        
        
        allCovs.save(opts.stats, r:true)
    }
    
    private Matrix loadCoverageStats(List<String> samples, List<String> covStatsFilePaths) {
        
        Matrix cov_stats =  Matrix.concat(
            samples.collect { sample ->
                Matrix.load(covStatsFilePaths.find { new File(it).name.startsWith(sample + '.') })
            }
        )
        cov_stats.sample = samples

        return cov_stats
    }

    private void writeCoeffVStats(String path, IntegerStats coeffVStats) {
        
        List cvThresholds = (0..100).step(5)
         
        List<List> percentiles = cvThresholds.collect { thresh ->
            [thresh, coeffVStats.fractionAbove(thresh)]
        }
       
        new File(path).withWriter { w ->
            w << 'coeffvPercentiles = // NOJSON\n'
            w << JsonOutput.prettyPrint(JsonOutput.toJson(percentiles))
            w << '\n'
        }
    }    

    private void writeCorrelationJS(String path, Matrix correlations) {
        List<String> all_samples = correlations.sample
        
        Map<String, Map<String,Double>> correlationMap = correlations.collect { row ->
           [
               sample,
               [ all_samples, row ].transpose().collectEntries()
           ]
        }.collectEntries()
        
        new File(path).withWriter { w ->
            w << 'corr = // NOJSON\n'
            w << JsonOutput.prettyPrint(JsonOutput.toJson(correlationMap))
            w << '\n'
        }
    }
    
    private void writeCovsJS(String path, Matrix cov_stats) {
        Map data = [
            means:  cov_stats.collect {
                [sample, delegate['Mean Coverage']]
            },
            medians: cov_stats.collect {
                [sample, delegate['Median Coverage']]
            },
        ]

        new File(path).withWriter { w ->
            w << 'corr = // NOJSON\n'
            w << JsonOutput.prettyPrint(JsonOutput.toJson(data))
            w << '\n'
        }
    }
    
    Matrix loadCovs(String path) {
        log.info "Load $path"
        def controls = Matrix.load(path)
        
        if(controls.properties.containsKey('Mean')) { // gngs bug
            controls.sample = controls.Mean
            controls.properties.remove('Mean')
        }
        
//        controls.names = (1..controls.columnDimension).collect { "C" + it}
//        controls.@displayColumns = 12
        return controls
    }
    
    static void main(String[] args) {
        cli('Calculates combined coverage statistics from multiple interval summary files', args) {
            corrTSV 'Output file to write row correlations to in TSV format', args:1, required: true
            corrJS 'Output file to write row correlations to in Javascript/JSON format', args:1, required: true
            covJS 'Output file to write combined coverage JS to', args:1, required: true
            stats 'Output file to write combined target coverage file to', args:1, required: true
            coeffvJS 'Output file to write the coefficient of variation percentiles to', args:1, required: true
            threads 'Number of threads to use in loading data', args:1, required: false, type: Integer
        }
    }
}
