
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
    
    /**
     * An (optional) set of true positives relevant to these results
     */
	CNVResults truth
}
