package ximmer.results
import groovy.lang.Closure;

import java.util.Map;

import gngs.RangedData
import gngs.Region


class DelfinResults extends CNVResults {
    

    DelfinResults(String fileName) {
        super(fileName, 0, 1, 2)
		this.load()
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        
        return super.load(options) { Region r ->
            r.sample = r.sample
            r.size = r.to - r.from
            r.quality = lr
            r.type = 'DEL'
        }
    }
}
