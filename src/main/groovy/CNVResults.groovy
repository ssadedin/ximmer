import groovy.json.JsonOutput

/**
 * Base class representing a file of CNV results
 * 
 * @author simon
 */
abstract class CNVResults extends RangedData {
    
	static Map RESULT_FACTORY = [
        'ed' : { new ExomeDepthResults(it) },
        'xhmm' : { new XHMMResults(it) },
        'cnmops' : { new CNMopsResults(it) },
        'ex' : { new ExcavatorResults(it) },
        'truth' : { new AngelResults(it) }
	]
    
	CNVResults(String fileName, int chrCol, int startCol, int endCol) {
        super(fileName, chrCol, startCol, endCol)
	}
    
	CNVResults(String fileName) {
        super(fileName)
	}
    
    String toJson(TargetedCNVAnnotator annotator = null) {
        JsonOutput.toJson(
                this.collect { cnv ->
                    def row = [
                        chr: cnv.chr,
                        start: cnv.from, 
                        end: cnv.to,
                        sample: cnv.sample,
                        quality: cnv.quality,
                        type: cnv.type
                    ]
                    
                    if(annotator) {
                        // Used to match the CNV type here, but after seeing that some
                        // CNV callers get confused and call a differentn type in a region where there is 
                        // a CNV,  I now think it's more accurate to return the frequency of ALL cnvs
                        // overlapping
                        row.spanningFreq = annotator.annotate(cnv, "ANY" /* cnv.type */).spanningFreq
                        row += annotator.annotateSize(cnv)
                    }
                    
                    if(truth != null) {
                        row.truth = truth.any { it.overlaps(cnv) && it.sample == cnv.sample }
                    }
                    
                    row
                }
        )
    }
    
    /**
     * An (optional) set of true positives relevant to these results
     */
	CNVResults truth
}
