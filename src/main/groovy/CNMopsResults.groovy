import groovy.lang.Closure;

import java.util.Map;

import gngs.Region


class CNMopsResults extends CNVResults {
    
    final static CN_TO_TYPE = [ 
            "CN0": "DEL", 
            "CN1" :"DEL",
            "CN2" :"NONE",
            "CN3" :"DUP",
            "CN4" :"DUP",
            "CN5" :"DUP",
            "CN6" :"DUP",
            "CN7" :"DUP"
            ]
    
    CNMopsResults(String fileName) {
        super(fileName)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        return super.load(options + [separator: ' ']) { Region r ->
            r.sample = r.sampleName
            r.quality = -r.median.toFloat()
            r.type = CN_TO_TYPE[r.CN]
        }
    }
    
    static void main(String [] args) {
       def cnmops = new CNMopsResults("testdata/11V00930_S1.merge.dedup.realign.recal.cn_mops_call_cnvs.tsv").load()
       
       println cnmops[0].sample + " region " + cnmops[0].toString() + " has " + cnmops[0].type
    }
}
