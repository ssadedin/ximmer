package ximmer

import gngs.ProgressCounter
import gngs.RangedData
import gngs.Region
import gngs.Regions
import gngs.ToolBase
import gngs.Utils
import graxxia.Matrix
import graxxia.Stats
import graxxia.TSV
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Use target read counts calculated in the QC stage of the Ximmer pipeline
 * to infer actual copy numbers for calls, along with reporting the observed
 * to expected read ratio and computing a combined quality metric that sums
 * the quality outputs from Phred based callers. These metrics are added
 * to Ximmer's standard TSV and JSON output formats.
 */
@Log
class InferCopyNumber extends ToolBase {
    
    static void main(String [] args) {
        cli('Infer copy number of merged called events in Ximmer pipeline', args) {
           c 'Target region coverage counts', args:1, required: true, type: File
           t 'TSV formatted report to update with copy number inferences', args:1
           j 'JSON formatted report to update with copy number inferferences', args:1
           to 'Output for TSV formatted report', args:1, type: File
           jo 'Output for JSON formatted report', args:1, type: File
        }
    }

    @Override
    public void run() {
        
        File covsFile = opts.c
        
        Matrix covs = Matrix.load(covsFile.path)
        covs.@displayColumns = 10
        covs.@displayPrecision = 2
        
        if(!opts.t && !opts.j)
            throw new IllegalArgumentException("Please provide at least one of -t or -j to specify CNVs to infer copy number for")
        
        log.info "Loaded coverage file with $covs.columnDimension targets and $covs.rowDimension samples from $covsFile"
        
        int region_index = 0
        Regions targets = covs.names.collect { new Region(it, index: region_index++) } as Regions
        
        log.info "Inferred $targets.numberOfRanges targets from coverage matrix"
        
        log.info "Normalising coverage ..."
        Matrix norm = covs.normaliseRows().normaliseColumns()

        log.info "Finished normalising coverage"
        
        // Load the CNVs
        Regions cnvs_j
        if(opts.j) {
            cnvs_j = new groovy.json.JsonSlurper().parseText(new File(opts.j).text).collect { 
                new Region(it.chr, it.start, it.end, * : it)
            }
            
            log.info "Loaded ${cnvs_j.numberOfRanges} CNV calls from JSON file"
        }
        
        Regions cnvs_t
        if(opts.t) {
            cnvs_t = new RangedData(opts.t).load()
            log.info "Loaded ${cnvs_t.numberOfRanges} CNV calls from TSV file"
        }
        
        // For now, just insist they have to be the same
        if(cnvs_j && cnvs_t && (cnvs_j.numberOfRanges != cnvs_t.numberOfRanges))
            throw new IllegalArgumentException("The size of the JSON format CNVs and TSV format CNVs do not match. Please ensure the same CNVs are provided in both files")
        
        List cnv_sources = []
        if(cnvs_j)
            cnv_sources.add(cnvs_j)
        if(cnvs_t)
            cnv_sources.add(cnvs_t)

        // Because they are the same size we can pair them
        List<List<Region>> paired = cnv_sources.transpose()
        
        ProgressCounter counter = new ProgressCounter(withRate: true, withTime:true, log:log)
        
        paired.each { List sources ->
            Region cnv = sources[0]

            // Find the targets overlapped by the caller
            def result = inferCopyNumber(norm, targets, cnv)

            sources*.coverage_ratio = result.coverageRatio
            sources*.combined_qual = result.combinedQuality
            sources*.copy_number = result.copyNumber

            counter.count()
        }
        counter.end()
        
        // Because of various manipulations, writing out the originally loaded data directly risks writing things
        // we didn't want to. So instead, read the files again and transfer the informationt to them
        // instead
        if(opts.to) {
            writeTSVOutput(cnvs_t)

        }
        
        if(opts.jo) {
            writeJSONOutput(cnvs_j)
        }
    }
    
    /**
     * Write back out CNVs provided in TSV format to specified file
     * 
     * @param cnvs_t
     */
    void writeTSVOutput(final Regions cnvs_t) {
        List<Map> rawTSV = new TSV(opts.t).toListMap()
        [rawTSV, cnvs_t].transpose().each { Map tsvRow, cnv ->
            tsvRow.putAll(
                coverage_ratio : cnv.coverage_ratio,
                combined_qual : cnv.combined_qual,
                copy_number : cnv.copy_number
            )
        }
        TSV.save(rawTSV, opts.to.path)        
        
        log.info "Wrote $opts.to.path"
    }
    
    /**
     * Write back out CNVs provided in JSON format, preserving line-oriented JSON form
     * of the original.
     * 
     * @param cnvs_j
     */
    void writeJSONOutput(Regions cnvs_j) {
        Utils.writer((File)opts.jo).withWriter { Writer w ->
            w.write('[\n')
            cnvs_j.eachWithIndex { cnv, i ->
                if(i > 0)
                    w.write(',\n')
                w.write(JsonOutput.toJson(cnv))
            }
            w.write('\n]')
        }        
        log.info "Wrote $opts.jo"
    }
    
    /**
     * Infer copy number and combined quality estimate for given CNV
     * 
     * @param cnvs
     */
    CNVCopyNumber inferCopyNumber(Matrix norm_covs, Regions indexedTargets, Region cnv) {

        def cnv_target_indices = indexedTargets.findIndexValues { it.overlaps(cnv) }
        int sample_index = norm_covs.sample.indexOf(cnv.sample)
        def cov_values = norm_covs[sample_index][cnv_target_indices]
        def stats = Stats.from(cov_values)
        
        List<String> callers = ['ed','xhmm','savvy']
        
        // Assume: all callers are in phred or similarly scaled values with a meaningful range of 0-100
        // so we are clipping values at each end of [0,100] range and then summing them
        double quality = callers.grep { cnv[it] == 'TRUE' }.collect {
            Math.min(100d, Math.max(0d, cnv[it + '_qual']))
        }.sum()
        
        // Currently copy number simply inferred by rounding to nearest integer closest based on 1.0 mean coverage ratio representing full diploid
        // quotient (2.0 copies). This will be wrong for males on the X chromosome
        int cn = Math.round(stats.mean * 2)
        
        CNVCopyNumber result = new CNVCopyNumber(
            coverageRatio : stats.mean,
            combinedQuality: quality,
            copyNumber : cn
        )
        return result
    }
}

@CompileStatic
class CNVCopyNumber {
    int copyNumber
    double coverageRatio
    double combinedQuality
}
