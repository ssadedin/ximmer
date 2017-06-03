/**
 * Utility functions
 */
function max(x,fn) {
  return x.reduce((maxVal, val) => {
    var txVal = (typeof(fn) === 'undefined') ? val : fn(val);
    return ((maxVal === null) || (txVal > maxVal)) ? txVal : maxVal;
  },null);
}

function min(x,fn) {
  return x.reduce((minVal, val) => {
    var txVal = (typeof(fn) === 'undefined') ? val : fn(val);
    return ((minVal === null) || (txVal < minVal)) ? txVal : minVal;
  },null);
}

function humanSize(value) {
    let units=['bp','kb','Mb','Gb']
    var fmt = d3.format('.1f');
    window.fmt = fmt;
    for(unit of units) {
        if(value < 1000) {
            console.log('value = ' + value + ' formatted = ' + fmt(value))
            return (fmt(value) ).replace(/.0$/,'') + unit
        }
        value = value / 1000
    }
}

/**
 * Half-open range class.
 * <p>
 * Range is inclusive of start, exclusive of end
 */
class Range {
    
    constructor(props, start, end) {
        // If only 1 args, assume object syntax
        if(!end) {
            Object.assign(this, props)
        }
        else {
            this.chr = props;
            this.from = start;
            this.to = end;
        }
    }
    
    containsWithinBounds(x) {
        return x >= this.from && x < this.to
    }
    
    overlaps(b) {
        let a = this;
        let result = (a.chr == b.chr) && 
                                 (a.containsWithinBounds(b.to) || 
                                  a.containsWithinBounds(b.from) ||
                                  b.containsWithinBounds(a.to))
        return result;
    }    
    
    toString() {
        return `${this.chr}:${this.from}-${this.to}`;
    }
}

/**
 * Genome chromosome sizes
 */

hg19_chr_sizes = {
    '1':249250621,
    '2':243199373,
    '3':198022430,
    '4':191154276,
    '5':180915260,
    '6':171115067,
    '7':159138663,
    'X':155270560,
    '8':146364022,
    '9':141213431,
    '10':135534747,
    '11':135006516,
    '12':133851895,
    '13':115169878,
    '14':107349540,
    '15':102531392,
    '16':90354753,
    '17':81195210,
    '18':78077248,
    '20':63025520,
    'Y':59373566,
    '19':59128983,
    '22':51304566,
    '21':48129895
}

/**
 * The frequency above which a CNV is not counted as a false positive due to
 * likely being a real CNV that is present in the population
 */
MAX_RARE_CNV_FREQ = 0.01

// ----------------------------- QScore Calibration ----------------------------------

/**
 * Ximmer CNV Evaluation Javascript
 * 
 * Displayes plots and summary information about CNV detection performance.
 */
class CallerCalibrationCurve {
  
  /**
   * calls - a list of objects with properties: id, caller, cnvs
   */
  constructor(calls) {
    this.calls = calls;
    this.width = 600,
    this.height = 260,
    this.margin = { bottom: 60, left: 60 };
  }
  
  calculateBins() {
    this.calls.forEach(caller => {
        caller.bins = this.calculateCallerBins(caller.cnvs);
        console.log(`Caller ${caller.id} has ${caller.bins.length} bins `)
    })
  }
  
  calculateCallerBins(cnvs) {
    
    const caller_max = max(cnvs, cnv => cnv.quality)
    const caller_min = min(cnvs, cnv => cnv.quality)

    let bin_size = (caller_max - caller_min) / 10
    bin_size = Math.round(bin_size / 5) * 5

    let bins = [];
    let binUpper = caller_min;
    let binMax = caller_max + bin_size;
    
    while(bin_size && (binUpper < binMax)) {
      bins.push({low: binUpper, high:binUpper+bin_size, count: 0, truth: 0})
      binUpper += bin_size
    }

    window.bins = bins;

    cnvs.forEach(cnv => { 
      var b = bins.find(b => cnv.quality > b.low && cnv.quality <= b.high); 
      if(b) {
        b.count++; 
       if(cnv.truth) b.truth++; 
      }
    });
    
    console.log(`Max quality = ${caller_max}, min quality = ${caller_min}, bin size = ${bin_size}, ${bins.length} bins`);

    // Remove bins that have fewer than 3 counts
    return bins.filter(bin => { return bin.count > 0 });
  }
  
  render() {
    throw "Please override the render method to provide a rendering implementation";
  }
}

class C3CallerCalibrationCurve extends CallerCalibrationCurve {
    
    constructor(calls) {
      super(calls);
    }
    
    render(id) {
      
       this.calculateBins();
      
       let xs = this.calls.reduce((result,callset) => {
         result[callset.caller + ' Precision'] = callset.id + '_x';
         return result;
       }, {});
      
       let xCols = this.calls.map(callset => {
                     return [callset.id + '_x'].concat(callset.bins.map(bin => (bin.low+bin.high)/2));
                   });
      
       let yCols = this.calls.map(callset => {
                     return [callset.caller + ' Precision'].concat(callset.bins.map(bin => bin.truth / bin.count));
                   });
           
       let cols = xCols.concat(yCols);
     
       c3.generate({
           bindto: '#' + id,
           data: {
               xs: xs,
               columns: xCols.concat(yCols),
               type: 'line'
           },
           grid: {
             x: {
               show: true
             },
             y: {
               show: true
             }
           },
           zoom: {
             enabled: true
           },
           axis: {
             y: { 
               label: {
                 text: 'Emprical Precision',
                 position: 'outer-middle'
               },
             },
             x: {
                 label: {
                     text: `${this.calls.map(c => c.caller).join(",")} Quality Scores`,
                     position: 'outer-center',
                     show:true,
                     style: 'font-size: 18px'
                 }
             }
         },
         legend: {
             position: 'right'
         }
       }); 
    }
}

function showQscores() {
    console.log("Showing qscore results");
    
    var callList = Object.keys(cnv_calls).map( caller => { return {
       id:  caller,
       'caller': caller,
       cnvs: cnv_calls[caller]
    }}).filter(calls => calls.id != "truth");
    
    let windowWidth = layout.center.state.innerWidth;
    let windowHeight = layout.center.state.innerHeight;
    let borderMargin = {x: 30, y:30};
    let plotMargin = {x: 30, y:30 }; // get from css somehow?
    
    // Size the plots into an even grid
    let minWidth = 300;
    let minHeight = 200;
    let maxHeight = 600;
    
    // Try to layout in 1 row at first
    let columns = callList.length;
    let calcMargin = () => plotMargin.x * columns + borderMargin.x*2; 
    let calcPlotWidth = () => Math.max(minWidth, Math.floor((windowWidth-calcMargin())/columns));
    
    // Wrap if exceeding the available width 
    let rows = 1;
    while(calcPlotWidth()*columns + calcMargin() > windowWidth) {
        let totalWidth =calcPlotWidth()*columns + calcMargin();
        console.log(`Columns = ${columns} rows = ${rows}, Total width: ${totalWidth} vs window: ${windowWidth}`);
        rows += 1;
        columns = Math.ceil(callList.length / rows);
        if(rows>5)
            break;
    }
    
    let plotWidth = calcPlotWidth();
    let plotHeight = Math.floor(Math.min(maxHeight,(windowHeight - plotMargin.y * rows - borderMargin.y) / rows));
    
    console.log(`Calculated ${columns}x${rows} grid for layout, plotWidth=${plotWidth}`)
    
    window.callList = callList;
    window.cc = new C3CallerCalibrationCurve([callList[0]])
    // cc.render('qscore_calibration_figure');
    
    with(DOMBuilder.dom) {
        callList.filter(calls => calls.id != 'truth').forEach(calls => {
            let plotId = 'ximmer_qscore_calibration_'+calls.id;
            
            let plot = DIV({ style:`display: inline-block; width: ${plotWidth}px; height: ${plotHeight}px;`},
                        H3({'class':'qualCalTitle'},'Quality Score Calibration Curve for ' + calls.id),
                        DIV({id: plotId}, calls.id)
                        );
            $('#qscore_calibration_figure')[0].appendChild(plot);
            new C3CallerCalibrationCurve([calls]).render(plotId);
        });
    }
}

// ----------------------------- Breakdown by Sample  ----------------------------------


class CNVsBySample {
    constructor(props) {
    }
    
    /**
     * Create a map indexed by caller, with values being a child map indexed by sample having 
     * values representing the total number of CNV calls for that sample.
     */
    calculateCounts(cnv_calls) {
        
        let countsBySample = Object.keys(cnv_calls).reduce((callerCounts,caller) => {
            if(caller == 'truth')
                return callerCounts;
            
            callerCounts[caller] = cnv_calls[caller].reduce((counts,cnv) => { 
                let sample = cnv.sample.replace(/-[^-]*$/,'');
                counts[sample] = counts[sample] ? counts[sample]+1 : 1; return counts; 
            }, {})
            return callerCounts;
        }, {})
        return countsBySample;
    }
    
    render(id) {
        
        let callers = Object.keys(cnv_calls).filter(c => c != 'truth');
        
        let counts = this.calculateCounts(cnv_calls);
        
        let samples = Object.values(counts).reduce((samples,callerCounts) => {
            return Object.keys(callerCounts).reduce((samples,sample) => { samples[sample] = true; return samples; }, samples)
        },{})
        
        // We actually only ever wanted a unique list of samples, so extract the keys
        // from the map
        samples = Object.keys(samples)
        
        let count_data = samples.map(function(sample) {
                return {
                    key: sample,
                    values: callers.map((caller, i) => { 
//                        return {label: caller, value: counts[caller][sample]}
                        return {series: caller, label: caller, x: i, y: counts[caller][sample]}
                    })
               }
        })
        
        window.count_data = count_data;
        
        var chart;
        nv.addGraph(function() {
            chart = nv.models.multiBarChart()
                .barColor(d3.scale.category20().range())
                .duration(300)
                .margin({bottom: 100, left: 70})
                .rotateLabels(45)
                .groupSpacing(0.1)
            ;

            chart.reduceXTicks(false).staggerLabels(true);

            chart.xAxis
                .axisLabel("Sample")
                .axisLabelDistance(35)
                .showMaxMin(false)
                .tickFormat((x) => callers[x])
            ;

            chart.yAxis
                .axisLabel("Count of CNV Calls")
                .axisLabelDistance(-5)
            ;

            d3.select('#'+id)
                .datum(count_data)
                .call(chart);

            nv.utils.windowResize(chart.update);

            return chart;
        });
    }
}


function showSampleCounts() {
    console.log("Showing sample count plot");
    let plot = new CNVsBySample(); 
    plot.render('cnvs_by_sample_chart')    
}


// ----------------------------- ROC Curve Functions  ----------------------------------

class CNVROCCurve {
    
    constructor(props) {
        /**
         * cnvs - a map keyed on caller id with values being the calls
         */
        this.rawCnvs = props.cnvs;
        this.maxFreq = props.maxFreq ? props.maxFreq : MAX_RARE_CNV_FREQ;
        this.sizeRange = props.sizeRange;
        
        if(!this.rawCnvs.truth) 
            throw new Error("ROC Curve requires true positives specified in CNV calls as 'truth' property")
    }
    
    computeROCStats(cnvs) {
        
        // the set of true positives that we have identified
        let tps = this.rawCnvs.truth.map(function(cnv, i) {
            var tp = Object.assign(new Range(cnv.chr, cnv.start, cnv.end),cnv);
            tp.id = i;
            tp.detected = false;
            return tp;
        });
        
        window.tps = tps;
        
        let tpCount = 0;
        let fpCount = 0;
        
        cnvs.forEach(function(cnv) {
            cnv.range = new Range(cnv.chr, cnv.start, cnv.end);
            if(cnv.truth) {
                // Which tp? we don't want to double count
                let tp = tps.find(tp => tp.overlaps(cnv.range) && (tp.sample == cnv.sample))
                
                if(!tp) {
                    console.log(`CNV marked as true but does not overlap truth set: ${cnv.chr}:${cnv.start}-${cnv.end}`);
                }
                else
                if(!tp.detected) {
                    tp.detected = true;
                    ++tpCount;
                }
            }
            else {
                ++fpCount;
            }
            cnv.tp = tpCount;
            cnv.fp = fpCount
        });
    }
    
    render(id) {
        
        let sizeMin = Math.pow(10, this.sizeRange[0]);
        let sizeMax = Math.pow(10, this.sizeRange[1]);
        
        const unfilteredCount = Object.values(this.rawCnvs).reduce((n,caller) => n+caller.length, 0);
        console.log(`There are ${unfilteredCount} raw cnv calls`);
        
        console.log("Filtering by spanningFreq < " + this.maxFreq + " size range = " + this.sizeRange);
        
        // First, filter by maxFreq since that makes everything else faster
        // then sort each cnv caller's CNVs in descending order of quality
        this.filteredCnvs = {};
        Object.keys(this.rawCnvs).filter(caller => caller != 'truth').forEach((caller) =>
            this.filteredCnvs[caller] = 
                this.rawCnvs[caller].filter(cnv => cnv.spanningFreq < this.maxFreq &&
                                                  (cnv.end - cnv.start > sizeMin) && 
                                                  (cnv.end - cnv.start < sizeMax) && 
                                                  ((simulationType != 'replace') || (cnv.chr == 'chrX' || cnv.chr == 'X')))
                                    .sort((cnv1,cnv2) => cnv2.quality - cnv1.quality)
        );
        
        let filteredTruth = 
            this.rawCnvs.truth.filter((cnv) => (cnv.end - cnv.start > sizeMin) && (cnv.end - cnv.start < sizeMax));
        
        let cnvCount = Object.values(this.filteredCnvs).reduce((n,caller) => n+caller.length, 0);
        console.log(`There are ${cnvCount} cnv calls after filtering by spanningFreq<${this.maxFreq}`);
        
        window.cnvs = this.filteredCnvs;
            
        // Now iterate through each caller's CNVs and compute the number of true and false positives
        Object.values(this.filteredCnvs).forEach((cnvs) => this.computeROCStats(cnvs));
        
        /*
        
        let xCols = Object.keys(this.filteredCnvs).map(caller => [caller+'_x'].concat(this.filteredCnvs[caller].map(cnv => cnv.fp)))
        
        let yCols = Object.keys(this.filteredCnvs).map(caller => [caller + '_y'].concat(this.filteredCnvs[caller].map(cnv => cnv.tp)))
        
        let xs = Object.keys(this.filteredCnvs).reduce(function(result, caller) { 
            result[caller+'_y'] = caller + '_x';
            return result
        }, {});
            
        c3.generate({
           bindto: '#' + id,
           data: {
               xs: xs,
               columns: xCols.concat(yCols),
               type: 'scatter'
           }, 
        });
        */
        
        
        let points = [];
        Object.keys(this.filteredCnvs).forEach(caller => points.push({
            values: this.filteredCnvs[caller].map(cnv => { return { x: cnv.fp, y: cnv.tp }}),
            key: caller
        }))
        
        window.points = points;
        
        var chart = nv.models.lineChart()
                             .margin({left: 100})  //Adjust chart margins to give the x-axis some breathing room.
                             .useInteractiveGuideline(true)  //We want nice looking tooltips and a guideline!
                             .showLegend(true)       //Show the legend, allowing users to turn on/off line series.
                             .showYAxis(true)
                             .showXAxis(true)
                             .padData(true)
                             .yDomain([0,filteredTruth.length])
                             .forceX([0])
                             
                        ;
                    
        chart.xAxis.axisLabel('False Positives')
        chart.yAxis.axisLabel('True Positives')
        
//        chart.tooltip.contentGenerator(function (obj) { console.log('called'); return JSON.stringify(obj)})        ;
//        var tooltip = chart.interactiveLayer;
//        tooltip.contentGenerator = function (d) { return "FUG"; };
//        chart.tooltip.contentGenerator= function(data, elem) {
//            elem.innerHTML = 'FOO';
//        };
//        
//        window.chart = chart;
        
        console.log("rendering to " + id);
        d3.select('#' + id)
          .datum(points)  
          .call(chart);  
        
        window.chart = chart;
//        nv.utils.windowResize(function() { chart.update() });        
    }
}

function showROCCurve() {
    console.log("ROC Curve");
    
    $('#cnv_roc_curve_container').html('<svg style="display: inline;" id=cnv_roc_curve></svg><div id=slider_label></div><div style="display: block; margin-left: 100px" id=slider></div>');
    
    let makePlot = (range) =>  {
        let plot = new CNVROCCurve({
            cnvs: cnv_calls,
            sizeRange: range
        });
        plot.render('cnv_roc_curve');
    };
        
    let labelFn = (range) => {
        $( "#slider_label" ).html("CNV Size Range: " + Math.pow(10,range[0]) + "bp - " + Math.pow(10,range[1])+"bp");
    };
    
    let initialRange = [0, 7];
    
    var redrawTimeout = null;
    
    $(function() {
        $("#slider").slider({
          range: true,
          values: initialRange,
          min: 0,
          max: 7,
          step: 1,
          slide: function( event, ui ) {
            console.log(ui.values);
            labelFn(ui.values);
            if(redrawTimeout)
                clearTimeout(redrawTimeout);
            redrawTimeout = setTimeout(() => {
                makePlot(ui.values);
              }, 2000);
          }
        }).width(450);
     });    
    
    labelFn(initialRange);
    makePlot(initialRange);
}

// ----------------------------- Genome Distribution  ----------------------------------

class CNVGenomeDistribution {
    constructor(props) {
        Object.assign(this, props)
        
        if(!this.bin_size) {
            this.bin_size = 10 * 1000 * 1000
            console.log("Using default bin size = " + this.bin_size)
        }
            
        if(!this.cnvCalls)
            this.cnvCalls = window.cnv_calls; // hack
    }
    
    
    computeBins() {
        // Create the chunks we want to scan. For a reasonable sized plot
        // 10mb produces about the right resolution
        this.chrStarts = []
        this.bins = Object.keys(hg19_chr_sizes).reduce((bins,chr) => {
            let pos = 0
            let max = hg19_chr_sizes[chr]
            
            this.chrStarts.push(bins.length)
            
            while(pos < max) {
                bins.push(new Range(chr, pos, pos + (this.bin_size-1)))
                pos += this.bin_size
            }
            return bins;
        }, [])
        
        return this.bins;
    }
    
    calculateCallerDistribution(bins,cnvs) {
        return bins.map(bin => {
            return cnvs.reduce((n, cnv) => bin.overlaps(cnv.range) ? n+1 : n, 0);
        })
    }
    
    /**
     * Create a map indexed by caller, with values being a child map indexed by sample having 
     * values representing the total number of CNV calls for that sample.
     */
    calculateDistribution(cnv_calls) {
        
        this.bins = this.computeBins()
        
        // Scan along the genome in 10mb chunks and find the number of CNV calls
        // for each caller
        let dist = {}
        Object.keys(cnv_calls).forEach(caller => {
            if(caller == 'truth')
                return
                
            dist[caller] =  this.calculateCallerDistribution(this.bins,cnv_calls[caller])
        })
        
        return dist
    }
    
    render(id, showChrs) {
        
        // If no specific chr defined, show all
        let chrs = hg19_chr_sizes
        if(showChrs) {
            chrs = {}
            showChrs.forEach(chr => chrs[chr] = hg19_chr_sizes[chr])
        }
        
        let caller_dists = this.calculateDistribution(this.cnvCalls)
        
        let points = []
        Object.keys(this.cnvCalls).forEach(caller => {
            if(caller == 'truth')
                return
                
            points.push({
                key: caller,
                values: caller_dists[caller].map((bin,i) => { return {x: i, y: bin, chr: this.bins[i].chr}})
                                            .filter(point => chrs[point.chr] )
            })
        })
        
        var chart = nv.models.lineChart()
                             .margin({left: 100})  //Adjust chart margins to give the x-axis some breathing room.
                             .useInteractiveGuideline(true)  //We want nice looking tooltips and a guideline!
                             .showLegend(true)       //Show the legend, allowing users to turn on/off line series.
                             .showYAxis(true)
                             .showXAxis(true)
                             .padData(true)
                             .forceY([0])
                             
                        ;
                    
        let xLabel = this.xLabel ? this.xLabel : 'Genome Position'
        chart.xAxis.axisLabel(xLabel + ' (' + humanSize(this.bin_size) + ' bins)')
                   .tickValues(this.chrStarts)
                   .tickFormat((i) => this.bins[i].from == 0 ? this.bins[i].chr : this.bins[i].chr + ':'+this.bins[i].from)
        
        chart.yAxis.axisLabel('Count of Overlapping CNVs')
        
        
        console.log("rendering to " + id);
        d3.select('#' + id)
          .datum(points)  
          .call(chart);  
        
        window.chart = chart;
        return chart;
    }
}


function showCNVGenomeDistribution() {
    
    let callsToShow = rare_calls;
    
    console.log("Showing genome distribution plot");
    let plot = new CNVGenomeDistribution({bin_size: 5 * 1000 * 1000, cnvCalls: callsToShow}); 
    let chart = plot.render('cnv_genome_dist')    
    
    chart.lines.dispatch.on('elementClick', function(info) {
        let pointInfo = info[0]
        let index = pointInfo.pointIndex;
        let chr = plot.bins[index].chr
        let subPlot = new CNVGenomeDistribution({bin_size: 500 * 1000, callsToShow, xLabel: 'Position in Chromosome ' + chr }); 
        subPlot.render('cnv_chr_dist',[chr])
    });
    
}

// ----------------------------- CNV Size Binning  ----------------------------------

class CNVSizePlot {
    constructor(props) {
        this.rawCnvs = props.cnvs;
    }
    
    render() {
    }
}
 
// ----------------------------- CNV Loading Functions ----------------------------------


function loadCnvs(callback, runsToLoad, results) {
    
    console.log("loadCnvs");
   
    var oldScriptElement = document.getElementById('cnvs_load_script');
    if(oldScriptElement)
        oldScriptElement.parentNode.removeChild(oldScriptElement);
    
    let run = runsToLoad.pop();
    
    console.log(`Loading cnvs from run ${run}`);

    const script = document.createElement("script");
    script.id = 'cnvs_load_script';
    script.src = run + '/' + analysisName + '/report/cnv_calls.js'
    script.async = true;
    
    let mergeResults =  () => {  
        Object.values(cnv_calls).forEach(cnvs => cnvs.forEach(cnv => { 
            cnv.sample = cnv.sample + "-" + run;
            cnv.range = new Range(cnv.chr.replace('chr',''), cnv.start, cnv.end);
        }));
        if(!results)
            results = cnv_calls;
        else
            Object.keys(results).forEach(caller => results[caller] = results[caller].concat(cnv_calls[caller]));
    };
    
    if(runsToLoad.length > 0) {
        script.onload = () => {
            mergeResults()
            loadCnvs(callback, runsToLoad, results);
        }; 
    }
    else {
        script.onload = () => {
            mergeResults()
            window.cnv_calls = results;
            window.rare_calls = {};
            
            // A number of different outputs are based on CNVs absent from the population,
            // so 
            Object.keys(cnv_calls).filter(caller => caller != 'truth').forEach((caller) => {
                window.rare_calls[caller] = cnv_calls[caller].filter(cnv => cnv.spanningFreq < MAX_RARE_CNV_FREQ)
            });
            callback()
        };
    }
    
    console.log("Loading cnv calls from " + script.src);
    document.body.appendChild(script); 
}

function loadCnvReport(runId, callback) {
    
    console.log("loadCnvReport");
   
    var oldScriptElement = document.getElementById('cnvs_load_script');
    if(oldScriptElement)
        oldScriptElement.parentNode.removeChild(oldScriptElement);
    
    const script = document.createElement("script");
    script.id = 'cnvs_load_script';
    script.src = runId + '/' + analysisName + '/report/cnv_report.b64.js'
    script.async = true;
    
    script.onload = callback
   
    console.log("Loaded cnv report from " + script.src);
    document.body.appendChild(script); 
}

function showCNVReport(runIndex, runId) {
    
    console.log('Showing CNV report for ' + runId + ' with ' + cnvReportHTML.length + ' bytes of HTML');
    
    var ifel = document.getElementById('run'+runIndex + 'Iframe');
    var doc = ifel.contentWindow.document;
    ifel.width = ($(window).width() - 50);
    ifel.height = ($(window).height() - 180);
    
    console.log("Write run " + runId + " to iframe");
    doc.open();
    doc.write(atob(cnvReportHTML));
    doc.close();                                
}

function loadAndCall(fn) {
    console.log("loadAndCall");
    
    if(typeof(window.cnv_calls) == 'undefined') {
        loadCnvs(fn, runs.map((r) => r)) // hack to clone runs
    }
    else {
        fn();
    }
}

$(document).ready(function() {
   console.log('Summary report init');
   $(events).on("activate", (event,ui) => {
       console.log("Activated in summary");
       var panelId = ui.newPanel[0].id;
       if(panelId == "qscorecalibration") {
           loadAndCall(showQscores);
       }
       else 
       if(panelId == "simrocs") {
           loadAndCall(showROCCurve);
       }
       else 
       if(panelId == "sample_counts") {
           loadAndCall(showSampleCounts);
       }
       else 
       if(panelId == "genome_dist") {
           loadAndCall(showCNVGenomeDistribution);
       } 
       else
       if(panelId.match(/runcalls[0-9]*/)) {
           let runIndex = panelId.match(/runcalls([0-9])*/)[1];
           
           let runId = runs[parseInt(runIndex,10)];
           
           console.log(`Show calls ${runId}`);
           
           loadCnvReport(runId, () => showCNVReport(runIndex, runId));
       }
   }); 
})