package ximmer

import gngs.*
import graxxia.*
import groovy.json.JsonOutput

import java.awt.Color
import java.util.logging.Logger
import gngs.plot.Bars
import gngs.plot.Plot
import gngs.plot.bx.*


/**
 * A tool that produces a basic report, mainly based on CNV counts to allow
 * quick assessment of whether the CNV calling has operated correctly.
 */
class CNVQCReport extends ToolBase {
    
    
    static Logger log = Logger.getLogger("CNVQCReport")

    Map<String,BigDecimal> chrSizes = [
        [ "chr": "chr1", "mb": 248.956422 ],
        [ "chr": "chr2", "mb": 242.193529 ],
        [ "chr": "chr3", "mb": 198.295559 ],
        [ "chr": "chr4", "mb": 190.214555 ],
        [ "chr": "chr5", "mb": 181.538259 ],
        [ "chr": "chr6", "mb": 170.805979 ],
        [ "chr": "chr7", "mb": 159.345973 ],
        [ "chr": "chr8", "mb": 145.138636 ],
        [ "chr": "chr9", "mb": 138.394717 ],
        [ "chr": "chr10", "mb": 133.797422 ],
        [ "chr": "chr11", "mb": 135.086622 ],
        [ "chr": "chr12", "mb": 133.275309 ],
        [ "chr": "chr13", "mb": 114.364328 ],
        [ "chr": "chr14", "mb": 107.043718 ],
        [ "chr": "chr15", "mb": 101.991189 ],
        [ "chr": "chr16", "mb": 90.338345 ],
        [ "chr": "chr17", "mb": 83.257441 ],
        [ "chr": "chr18", "mb": 80.373285 ],
        [ "chr": "chr19", "mb": 58.617616 ],
        [ "chr": "chr20", "mb": 64.444167 ],
        [ "chr": "chr21", "mb": 46.709983 ],
        [ "chr": "chr22", "mb": 50.818468 ],
        [ "chr": "chrX", "mb": 156.040895 ],
        [ "chr": "chrY", "mb": 57.227415 ]
    ].collectEntries { [it.chr, it.mb] }

    static void main(String [] args) {
        cli('SVQCReport -sample <sample id> <cnvs json file>', args) {
            sample 'Name of sample', args: 1, required: true
            batch 'Name of batch sample was run in', args: 1, required: false
            json 'Output JSON file', args: 1, required: true
            o 'Output PDF file name', args:1, required: true
        }
    }

    @Override
    public void run() {

        String sample = opts.sample

        Regions cnvs = RangedData.loadJSON(opts.arguments()[0])

        log.info "Loaded ${cnvs.numberOfRanges} CNVs"
        
        // Take advantage of fact that callers are each listed in the keys of each row
        // with their quality scores to infer the set of callers
        List callers = cnvs[0].properties*.key.grep { it.endsWith('_qual') }.collect { it.tokenize('_')[0] }
        
        log.info "Inferred callers as " + callers

        
        Map qcJSON = [
            Sample: sample,
            counts: [
                total: cnvs.numberOfRanges
            ],
            callers : [:]
        ]
        
        if(opts.batch)
            qcJSON.batch = opts.batch

        log.info "Creating CNVS size distribution plot"
        Plot sizePlot = new Plot(title: "CNV Size Distribution for $opts.sample", xLabel: 'CNV Size', yLabel: 'Frequency') << \
            new Density.Area(data: cnvs.grep { it.size() < 20000 }*.size())
            
        sizePlot.save("${sample}_cnv_size_distribution.png")

        log.info "Creating CNVS small size distribution plot"
        sizePlot = new Plot(title: "Small CNV Size Distribution for $opts.sample", xLabel: 'CNV Size', yLabel: 'Frequency') << \
            new Density.Area(data: cnvs.grep { it.size() < 2000 }*.size())
        sizePlot.save("${sample}_small_cnv_size_distribution.png")


        log.info "Creating Per Caller size distribution plot"
        sizePlot = new Plot(title: "Per-Caller CNV Size Distribution for $opts.sample", xLabel: 'CNV Size', yLabel: 'Frequency')
        for(caller in callers) {
            sizePlot << new Density.Area(data: cnvs.grep { it[caller] == "TRUE" && it.size() < 5000 }*.size(), displayName: caller)
        } 
        sizePlot.save("${sample}_per_caller_cnv_size_distribution.png")


        Map<String,Integer> chrCounts = cnvs.countBy { it.chr }
        qcJSON.counts.by_chr = chrCounts

        Plot chrPlot = new Plot(title: "Count of CNVs by Chromosome", xLabel: "Chromosome", yLabel: "Count") << \
            new Bars(x: (1..chrCounts.size()), y: chrCounts*.value, labels: chrCounts*.key)
            
        chrPlot.save("${sample}_cnv_counts_by_chr.png")

        chrPlot = new Plot(title: "Count of CNVs by Chromosome / Mb", xLabel: "Chromosome", yLabel: "Count") << \
            new Bars(x: (1..chrCounts.size()), y: chrCounts.collect { it.value  / chrSizes[it.key]}, labels: chrCounts*.key)
            
        chrPlot.save("${sample}_cnv_counts_by_chr_per_mb.png")

        log.info "Creating PDF $opts.o"
        PDF pdf = new PDF().document(opts.o) {
            
            font(Color.black, 16) {
                bold {
                    p("CNV QC Report for $sample")
                    if(opts.batch)
                        p("Batch: $opts.batch")
                }
            }

            table {
                head {
                    cell("Caller")
                    cell("Total Count")
                    cell("Solo Count") 
                    cell("Solo %") 
                }

                callers.each { caller ->
                    int totalCount = cnvs.count {  it[caller] == "TRUE" }
                    int soloCount = cnvs.count {  it[caller] == "TRUE" && it.count == 1 }
                    cell(caller) 
                    cell(totalCount)
                    cell(soloCount)
                    
                    double soloFrac = soloCount / (double)totalCount
                    qcJSON.callers[caller] = [
                        total: totalCount,
                        solo: soloCount,
                        soloFrac: soloFrac
                    ]
                    
                    color(soloFrac < 0.4 ? "black" : "red") {
                        cell(Utils.perc(soloFrac))
                    }
                }
            }
            
            img("${sample}_cnv_size_distribution.png")

            img("${sample}_small_cnv_size_distribution.png")

            img("${sample}_per_caller_cnv_size_distribution.png")
            
            img("${sample}_cnv_counts_by_chr.png")
            
            img("${sample}_cnv_counts_by_chr_per_mb.png")
        }
        
        new File(opts.json).withWriter { w ->
            w << JsonOutput.prettyPrint(JsonOutput.toJson(qcJSON))
            
        }
        
        log.info "Wrote $opts.json"
        
        log.info "Done"
    }
}

