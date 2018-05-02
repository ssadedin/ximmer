import groovy.lang.Closure;

import java.util.Map;

import gngs.RangedData


class PositiveControlResults extends CNVResults {
    
    PositiveControlResults(String fileName) {
        super(fileName)
    }

    @Override
    public RangedData load(Map options, Closure c = null) {
        return super.load(options + [columnNames: ['chr','start','end','sample']]) { r ->
            r.type = "DEL"
            r.quality = 1000
        }
    }
}
