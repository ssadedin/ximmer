/**
 * A set of VueJS components to display QC plots for CNV data
 */

// ----------------------------- QScore Calibration ----------------------------------

define(['vue','nv.d3'], function(Vue, nvd3) {

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
        this.minBinCount = 3
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
        
        console.log("caller max = " + caller_max + " caller min " + caller_min)
        
        let bin_size = (caller_max - caller_min) / 5
        bin_size = (Math.round(10 * bin_size / 5) * 5) / 10;
    
        let bins = [];
        let binUpper = caller_min;
        let binMax = caller_max + bin_size;
        
        while(bin_size && (binUpper < binMax)) {
          bins.push({low: binUpper, high:binUpper+bin_size, count: 0, truth: 0})
          binUpper += bin_size
        }
    
        window.bins = bins;
    
        cnvs.forEach(cnv => { 
            
          // Skip CNVs that are in DGV
          if(cnv.DGVFreq > MAX_RARE_CNV_FREQ) 
              return
              
          if((simulationType == 'replace') && (cnv.chr != 'X') && (cnv.chr != 'chrX'))
              return
            
          var b = bins.find(b => cnv.quality >= b.low && cnv.quality < b.high); 
          
          if(b) {
            b.count++; 
           if(cnv.truth) b.truth++; 
          }
        });
        
        console.log(`Max quality = ${caller_max}, min quality = ${caller_min}, bin size = ${bin_size}, ${bins.length} bins`);
    
        // Remove bins that have fewer than 3 counts
        let candidateBins =  bins.filter(bin => { return bin.count >= this.minBinCount });
        
        // if we ended up with only 1 bin, take the CNVs from that bin and re-bin
        if(candidateBins.length == 1) {
            let bin = candidateBins[0];
            console.log("Naive binning produced too few bins: exploding center bin from ${bin.low}-${bin.high}")
            let newBins = this.calculateCallerBins(cnvs.filter(cnv => cnv.quality >= bin.low && cnv.quality <= bin.high ))
            return newBins;
        }
        
        return candidateBins;
      }
      
      render() {
        throw "Please override the render method to provide a rendering implementation";
      }
    }
    
    class NVD3CallerCalibrationCurve extends CallerCalibrationCurve {
        constructor(calls) {
          super(calls);
        }
        
       render(id) {
           this.calculateBins();
           
           let values = this.calls[0].bins.map(bin => {
               return { 
                   x: (bin.low + bin.high)/2,
                   y: bin.truth/bin.count
               }
           })
           
           let points = [{
               key: this.calls[0].caller + " Quality Scores",
               values: values
           }];
           
            let binEdges = this.calls[0].bins.map(bin => Math.round(bin.low)).concat([Math.ceil(this.calls[0].bins[this.calls[0].bins.length-1].high)])
            
            var chart = nv.models.lineChart()
                                 .margin({left: 100})  //Adjust chart margins to give the x-axis some breathing room.
                                 .useInteractiveGuideline(true)  //We want nice looking tooltips and a guideline!
                                 .showLegend(true)       //Show the legend, allowing users to turn on/off line series.
                                 .showYAxis(true)
                                 .yDomain([0,1])
                                 .forceY([0,1])
                                 .showXAxis(true)
                                 .interpolate('basis')
                                 .padData(true)
                                 .forceX(binEdges)
                                 
                            ;
                        
            
            
            console.log("Bin edges are: " + binEdges.join(','))
            
            chart.xAxis.axisLabel('Quality Scores')
                       .tickValues(binEdges)
                       
            chart.yAxis.axisLabel('Empirical Precision')
                       .tickValues([0,0.2,0.4,0.6,0.8,1.0, 1.2])
            
            d3.select('#'+id)
                .datum(points)
                .call(chart);
            
            return chart;
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
        let maxHeight = 500;
        
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
        
        let heightFactor = 2.0;
        
        let plotWidth = calcPlotWidth();
        let plotHeight = Math.floor(Math.min(maxHeight,(windowHeight - plotMargin.y * rows - borderMargin.y) / rows));
        
        console.log(`Calculated ${columns}x${rows} grid for layout, plotWidth=${plotWidth}`)
        
        window.callList = callList;
    //    window.cc = new C3CallerCalibrationCurve([callList[0]])
        // cc.render('qscore_calibration_figure');
        
        with(DOMBuilder.dom) {
            callList.filter(calls => calls.id != 'truth').forEach(calls => {
                let plotId = 'ximmer_qscore_calibration_'+calls.id;
                let plotIdWrapper = plotId + '_wrapper'
                console.log("Append " + plotIdWrapper)
                
                let plot = DIV({ style:`display: inline-block; width: ${plotWidth}px; height: ${plotHeight}px;`},
                            DIV({id: plotIdWrapper})
                            );
                
                $('#qscore_calibration_figure')[0].appendChild(plot);
                
                $('#' + plotIdWrapper).html('<svg id="' + plotId + '" + style="' + `display: inline-block; width: ${plotWidth}px; height: ${plotHeight}px;`+'"></svg>')
                
                new NVD3CallerCalibrationCurve([calls]).render(plotId);
            });
        }
    }
    
    // ----------------------------- Breakdown by Sample  ----------------------------------
    
    
    class CNVsBySample {
        constructor(cnv_calls) {
            this.cnv_calls = cnv_calls
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
            
            let callers = Object.keys(this.cnv_calls).filter(c => c != 'truth');
            
            let counts = this.calculateCounts(this.cnv_calls);
            
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
            var self = this;
    //        nv.addGraph(function() {
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
    
                chart.multibar.dispatch.on('elementClick', (info) => {
                    let callerIndex = info.index;
                    let callerId = callers[callerIndex]
                    
                    let sampleIndex = info.data.series;
                    let sampleId = count_data[sampleIndex].key
                    
                    console.log('Hiding sample ' + sampleId + ' for caller ' + callerId)
                    
                    let samples_with_runs = runs.map(run_id => sampleId + '-' + run_id)
                    
                    let exclusions = model.excludedSamples[callerId] || []
                    if(exclusions.indexOf(samples_with_runs[0])<0) {
                        for(var s of samples_with_runs)
                            exclusions.push(s)
                    }
                    else {
                        exclusions = exclusions.filter(sampleToToggle => samples_with_runs.indexOf(sampleToToggle)<0)
                    }
                    model.excludedSamples[callerId] = exclusions
                    model.excludedSamples = Object.assign({}, model.excludedSamples)
                    self.cnv_calls = model.filterCNVs()
                    setTimeout(() => self.render(id), 300)
                }); 
                
                nv.utils.windowResize(chart.update);
    
                return chart;
    //        });
        }
    }

    Vue.component('sample-counts', {
            
        props: ['model'],
            
        methods: {
            showSampleCounts: function() {
                console.log("Showing sample count plot");
                let plot = new CNVsBySample(model.filterCNVs()); 
                plot.render('cnvs_by_sample_chart')    
            }
        },
            
        data: function() { return this.model },
            
        template: `
           <div>
            <h2>By Sample</h2>
            <div id=cnvs_by_sample_chart_container>
                <svg id=cnvs_by_sample_chart>
                </svg>
            </div>
          </div> 
        `
    })
        
    Vue.component('genome-distribution',{
            
        props: ['model'],
            
        data: function() { return this.model },
            
        methods: {
            showCNVGenomeDistribution: function() {
                
                let callsToShow = model.filterCNVs();
                    
                console.log("Showing genome distribution plot (vue)");
                let plot = new CNVGenomeDistribution({bin_size: 5 * 1000 * 1000, cnvCalls: callsToShow}); 
                let chart = plot.render('cnv_genome_dist')    
                    
                if(this.subPlot) {
                    console.log("Showing genome dist subplot (vue)");
                    this.subPlot.cnvCalls = callsToShow
                    this.subPlot.render('cnv_chr_dist',[this.subPlotChr])
                }
                    
                chart.lines.dispatch.on('elementClick', (info) => {
                    let pointInfo = info[0]
                    let index = pointInfo.pointIndex;
                    this.subPlotChr = plot.bins[index].chr
                    this.subPlot = new CNVGenomeDistribution({bin_size: 500 * 1000, cnvCalls: callsToShow, xLabel: 'Position in Chromosome ' + this.subPlotChr }); 
                    this.subPlot.render('cnv_chr_dist',[this.subPlotChr])
                });
            }
        },
            
        template: `
          <div>
            <h2>Distribution Along Genome</h2>
            <p>It is common that certain regions of the genome can present extreme difficulty in CNV calling and result in
               clusters of false positives that undermine the results. Click anywhere in the main plot to zoom to the specific
               chromosome in a sub-plot.
            </p>
            <div id=cnv_genome_dist_container>
                <svg id=cnv_genome_dist>
                </svg>
                <svg id=cnv_chr_dist>
                </svg> 
            </div>
          </div>
        `
    })
        
    Vue.component('sens-by-size', {
            
        props: ['model'],
            
        data: function() { return this.model },
          
        template: `
            <div>
                <h2>Performance by Deletion Size</h2>
                                
                <svg id=cnv_size_breakdown>
                </svg>
            </div>
        `
    })
        
    Vue.component('qual-score-calibration', {
            
        props: ['model'],
            
        data: function() { return this.model },
          
        template: `
            <div>
                <h2>Quality Score Calibration Curves</h2>
                <p>
                    These curves estimate the precision as quality varies for each CNV caller. Use these curves
                    to determine the appropriate level of filtering to control your false discovery rate.
                </p>
                <div id='qscore_calibration_figure'>
                </div>
            </div>
        `
    })

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
    
    // ----------------------------- CNV Size Binning  ----------------------------------
    
    class CNVSizePlot {
        constructor(props) {
            
            this.rawCnvs = props.cnvs;
            this.sizeParam = props.sizeParam; // eg: targets
            
            // Since different types of params can be specified, 
            // it is hard to create good bins automatically.
            // So for this case we just let the user specify them
            this.paramBins = props.paramBins;
            this.paramDesc = props.paramDesc;
            this.xScale = props.xScale  || function(x) { return x; };
            this.tickFormat = props.tickFormat
        }
        
        initBins() {
            
            let binCounts = []
            
            this.paramBins.map((bin,i) => {
                if(i > 0) {
                    binCounts.push({
                        index: i-1,
                        min: this.paramBins[i-1],
                        max: bin,
                        tp: 0,
                        cnvs: []
                    })
                }
            })
            
            return binCounts;
        }
        
        calculateBinsForCaller(cnvs) {
            let bins = this.initBins()
            cnvs.forEach( cnv => {
                // Add it to the first bin where it fits
                for(let i=0; i<bins.length; ++i) {
                    
                    let size = cnv[this.sizeParam]
                    if((size>=bins[i].min) && (size<bins[i].max)) {
                        bins[i].cnvs.push(cnv);
                        if(cnv.truth)
                                ++bins[i].tp;
                        break
                    }
                }
            })
            return bins;
        }
        
        /**
         * For a given bin containing true positive CNVs and a specified
         * CNV caller, compute the sensitivity for true positive CNVs in that bin.
         */
        calculateBinSensitivityForCaller(bin, caller) {
            // Copy the tps
            let tps = bin.cnvs.map(x => { return { cnv: x, detected: false}})
            if(tps.length == 0) {
                return 0
            }
                    
            let countTp = 0
                    
            let callerCnvs = cnv_calls[caller]
                    
            callerCnvs.forEach(cnv => {
                        
                if(!cnv.truth)
                    return
                        
                // Which tp? Might already have been detected
                let tp = tps.find(tp => tp.cnv.range.overlaps(cnv.range) && (tp.cnv.sample == cnv.sample))
                if(!tp) {
                    return
                }
                        
                if(!tp.detected) {
                    ++countTp
                    tp.detected = true
                }
            })
            
            return countTp / tps.length;
        }
        
        calculateSensititivyByBin(callerBins, callers) {
            callerBins.truth.forEach((bin, i) => {
                bin.callerSens = {}
                callers.forEach( caller => {
                    bin.callerSens[caller] = this.calculateBinSensitivityForCaller(bin,caller)
                })
            })
        }
        
        render(id) {
            
            let callers = Object.keys(cnv_calls)
            
            // find the number of true positives by different size properties
            let callerBins = {};
            callers.forEach(caller => callerBins[caller] = this.calculateBinsForCaller(cnv_calls[caller]));
            
            let realCallers = callers.filter(c => c != 'truth')
            
            this.calculateSensititivyByBin(callerBins, realCallers)
            
            let filteredBins = callerBins.truth.filter(bin => bin.cnvs.length>2)
            
            let minTruthCNVs = 3 
            
            let points = realCallers.map( (caller,callerIndex) => { 
                
                let values = [{x:0, y:0}].concat(filteredBins.map((bin,i) => { 
                                 return { 
                                     x: this.xScale((bin.min + bin.max)/2),
                                     y: bin.callerSens[caller]
                                 }
                            }))
                
                return { 
                    key: caller, 
                    values: values  
                }
            })
            
            let xMax = max(filteredBins.map(b => b.max))
            
            // Each bin in the binned truth set now has a callerSens property which is a 
            // map of caller => sensitiivty
            
            let chart = nv.models.lineChart()
                                 .margin({left: 100})  //Adjust chart margins to give the x-axis some breathing room.
                                 .useInteractiveGuideline(true)  //We want nice looking tooltips and a guideline!
                                 .showLegend(true)       //Show the legend, allowing users to turn on/off line series.
                                 .yDomain([0,1.05])
                                 .showYAxis(true)
                                 .showXAxis(true)
                                 .forceY([0])
                                 .forceX([this.xScale(xMax)])
                                 .pointShape('circle')
                                 .interpolate('basis')
                    ;
            
            
    //        if(this.scale)
    //            chart.xScale(this.scale);
      
            let xAxisLabel = this.paramDesc || this.sizeParam;
            
            chart.xAxis.axisLabel(xAxisLabel)
                       .tickValues(this.chrStarts)
                       
            if(this.tickFormat)
                   chart.xAxis.tickFormat(this.tickFormat)
              
            chart.yAxis.axisLabel('Fraction of True Positives Detected')
                       .tickValues([0,0.2,0.4,0.6,0.8,1.0])
            
            d3.select('#' + id)
              .datum(points)  
              .call(chart);  
          
            window.chart = chart;        
            return chart;
        }
    }
    
    function showSizeBreakdown() {
        
        let sizeTypes = [
            {
                sizeParam: 'targets',
                paramDesc: 'Number of Target Regions',
                paramBins: [0,1,2,3,4,5,10,20,50,100,500,1000]
            },
            {
                sizeParam: 'targetBp',
                paramDesc: 'Number of targeted Base Pairs',
                paramBins: [0,100,500,1000,5000,10000,50000,100000,500000,1000000,10000000],
            },
            {
                sizeParam: 'size',
                paramDesc: 'Genomic Span of CNV (bp)',
                paramBins: [0,100,500,1000,5000,10000,50000,100000,500000,1000000,10000000],
                xScale: Math.log10,
                tickFormat: (x) => '10^' + x
            }, 
        ]
        
        let showPlot = function(params) {
            let sizeChart = new CNVSizePlot(Object.assign({
                cnvs: cnv_calls
            }, params))
            sizeChart.render('cnv_size_breakdown');
        }
        
        if(!$('#size_radios').length) {
            with(DOMBuilder.dom) {
                let radios = DIV({id:'size_radios'},sizeTypes.map((st,index) => {
                    return DIV({},INPUT({type:'radio', value: index, name:'sizeradio'}),SPAN(st.paramDesc))
                }));
                
                $('#sizebreakdown')[0].appendChild(radios);
                $('#size_radios input')[0].checked = true;
                $('#size_radios input').change(function() {
                   console.log("Showing param set " + this.value);
                   showPlot(sizeTypes[parseInt(this.value,10)])
                });
            }
        }
        
        $('#size_radios input')[0].checked = true;   
        let sizeChart = new CNVSizePlot(sizeTypes[0])
        
        sizeChart.render('cnv_size_breakdown');
    }
    
    console.log("Loaded QC components")
    
    return {
        showSizeBreakdown : showSizeBreakdown,
        showQscores : showQscores
    }
})