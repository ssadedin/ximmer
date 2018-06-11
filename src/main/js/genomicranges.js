/*
 * Genomic Regions Indexing 
 * 
 * Shamelessley adapted from:
 * 
 * https://github.com/ucscXena/static-interval-tree
 */

/**
 * A set of ranges all belonging to a single contig of a genome
 */
class Ranges {
    
    constructor(all) {
        this.rangeIndex = this.index(all)
    }

    // Build a balanced binary tree from an ordered array.
     toTree(arr, low, high) {
        if (low >= high) {
            return undefined;
        }
        var mid = Math.floor((high + low) / 2);
        return {
            el: arr[mid],
            right: this.toTree(arr, mid + 1, high),
            left: this.toTree(arr, low, mid)
        };
    }
    
    getHigh({high}={high: -Infinity}) {
        return high
    }
    
    // Find the highest end value of each node. Mutates its input.
    findEnd(node) {
        if (!node) {
            return undefined;
        }
        var {left, right, el} = node;
        this.findEnd(left);
        this.findEnd(right);
        node.high = Math.max(this.getHigh(left), this.getHigh(right), el.end);
        return node;
    }
    
    // Build index.
    // intervals :: [{start, end, ...}, ...]
    index(intervals) {
        let cmp = (x, y) => x === y ? 0 : (x < y ? -1 : 1);
        var sorted = intervals.slice(0).sort((a, b) => cmp(a.start, b.start));
        return this.findEnd(this.toTree(sorted, 0, sorted.length));
    }
    
    matchesAcc(node, pos, acc) {
        if (node) {
            var {start, end} = pos;
            if (node.high >= start) {
                if (end >= node.el.start) {
                    if (node.el.end >= start) {
                        acc.push(node.el);
                    }
                    this.matchesAcc(node.right, pos, acc);
                }
                this.matchesAcc(node.left, pos, acc);
            }
        }
        return acc;
    }
    
    // Find intervals in node overlapping position.
    // pos :: {start, end}
    matches(pos) {
        return this.matchesAcc(this.rangeIndex, pos, [])
    }
    
    matches01Acc(node, pos, acc) {
        if (node) {
            var {start, end} = pos;
            if (node.high > start) {
                if (end > node.el.start) {
                    if (node.el.end > start) {
                        acc.push(node.el);
                    }
                    this.matches01Acc(node.right, pos, acc);
                }
                this.matches01Acc(node.left, pos, acc);
            }
        }
        return acc;
    }
    
    // Find intervals in node overlapping position, using half-open coords.
    // We could also support this by parameterizing the compare fn, though that adds
    // more function calls to what is expected to be a hot loop.
    // pos :: {start, end}
    matches01(node, pos) {
        return matches01Acc(node, pos, []);
    }
}

/**
 * A set of genomic regions, possibly belonging to heterogeneous contigs
 */
class Regions {
    constructor(all) {
        this.rangeLists = all.reduce((acc, range) => {
            acc[range.chr] = acc[range.chr] || [];
            acc[range.chr].push(range)
            return acc
        },{})
        
        this.chrs = Object.keys(this.rangeLists)
        this.ranges = new Map(this.chrs.map(chr => {
            return [chr, new Ranges(this.rangeLists[chr])]
        }))
    }
    
    getOverlaps(region) {
        let chrRanges = this.ranges.get(region.chr)
        return chrRanges ? chrRanges.matches(region) : []
    }
}

if(typeof process != 'undefined') {
    
    console.log('Running Tests')
    
    let out = console.log
    
    let r = new Regions([
        {chr:'chr1', start: 100, end: 200},
        {chr:'chr2', start: 120, end: 220},
        {chr:'chr1', start: 100, end: 200},
        {chr:'chr2', start: 150, end: 180},
    ])
    
    let o = r.getOverlaps({chr:'chr1',start: 102, end: 108})
    if(o.length != 2) 
            throw "Wrong overlaps for chr1:102-108: " + o.toString() 
            
    o = r.getOverlaps({chr: 'chr2', start: 124, end: 128})
    if(o.length != 1) 
            throw "Wrong overlaps for chr2:124-128: " + o.toString()                 
            
    out('Tests passed')
}

        