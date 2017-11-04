
// ----------------------------- ROC Curve Functions  ----------------------------------

class CNVROCCurve {
    
    constructor(props) {
        /**
         * cnvs - a map keyed on caller id with values being the calls
         */
        this.rawCnvs = props.cnvs;
        this.maxFreq = props.maxFreq ? props.maxFreq : MAX_RARE_CNV_FREQ;
        this.sizeRange = props.sizeRange;
        this.targetRange = props.targetRange;
        this.callSizeRange = props.callSizeRange;
        this.callTargetRange = props.callTargetRange;
        
        if(!this.targetRange) 
            this.targetRange = [0, 1000000];
        
        if((this.targetRange[1] == "Infinity") || (this.targetRange[1]<0)) {
            console.log("Infinite no. targets")
            this.targetRange[1] = 1000000;
        }
        
        if(!this.rawCnvs.truth) 
            throw new Error("ROC Curve requires true positives specified in CNV calls as 'truth' property")
    }
    
    computeROCStats(cnvs) {
        
        // the set of true positives that we have identified
        let tps = this.rawCnvs.truth.map(function(cnv, i) {
            var tp = {
                range: cnv.range,
                sample: cnv.sample,
                id: i,
                detected: false
            }
            return tp;
        });
        
        window.tps = tps;
        
        let tpCount = 0;
        let fpCount = 0;
        
        cnvs.forEach(function(cnv) {
            if(cnv.truth) {
                // Which tp? we don't want to double count
                let tp = tps.find(tp => tp.range.overlaps(cnv.range) && (tp.sample == cnv.sample))
                
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
                if(cnv.spanningFreq < MAX_RARE_CNV_FREQ)
                    ++fpCount;
            }
            cnv.tp = tpCount;
            cnv.fp = fpCount
        });
    }
    
    computeCallerLabelMap() {
        // Ximmer replaces the original configs with short versions ... but here we have to convert them
        // back to look up the labels
        let convertedConfigIds = analysisConfig.callerCfgs.map(cfg => { 
            let longId = cfg.split('_')[0]; 
            return cfg.replace(new RegExp('^' + longId), callerIdMap[longId])
        })
        
        return new Map(_.zip(convertedConfigIds, analysisConfig.callerLabels))
    }
    
    render(id) {
        
        this.filteredCnvs = model.filterCNVs()
        
        let rawCnvs = this.rawCnvs;
        let sizeMin = Math.pow(10, this.sizeRange[0]);
        let sizeMax = Math.pow(10, this.sizeRange[1]);
        let targetsMin = this.targetRange[0];
        let targetsMax = this.targetRange[1];
        
        let filteredTruth = 
            this.rawCnvs.truth.filter((cnv) => (cnv.targets >= targetsMin) && (cnv.targets<=targetsMax) && 
                                               (cnv.end - cnv.start > sizeMin) && (cnv.end - cnv.start < sizeMax) &&
                                               ((simulationType == 'replace') || (cnv.chr != 'chrX' && cnv.chr != 'X')))
                                               ;
        
        let cnvCount = Object.values(this.filteredCnvs).reduce((n,caller) => n+caller.length, 0);
        console.log(`There are ${cnvCount} cnv calls after filtering by spanningFreq<${this.maxFreq}`);
        
        window.cnvs = this.filteredCnvs;
            
        // Now iterate through each caller's CNVs and compute the number of true and false positives
        Object.values(this.filteredCnvs).forEach((cnvs) => this.computeROCStats(cnvs));
        
        let callerLabels = this.computeCallerLabelMap()
       
        let points = [];
        Object.keys(this.filteredCnvs).forEach(caller => points.push({
            values: this.filteredCnvs[caller].map(cnv => { return { x: cnv.fp, y: cnv.tp, quality: cnv.quality }}),
            key: callerLabels.get(caller)
        }))
        
        window.points = points;
        
        var chart = nv.models.lineChart()
                             .margin({left: 100, right: 130})  // note: right margin is mainly just to allow for the word 'Sensitivity'
                             .showLegend(true)       //Show the legend, allowing users to turn on/off line series.
                             .showYAxis(true)
                             .showXAxis(true)
                             .padData(true)
                             .yDomain([0,filteredTruth.length])
                             .forceX([0])
                        ;
                    
        chart.xAxis.axisLabel('False Positives')
        chart.yAxis.axisLabel('True Positives')
        
        let percFormat = d3.format('%0.1f')
        let fracFormat = d3.format('0.1f')
        
        console.log("Filtered truth has " + filteredTruth.length + " true positives")
        
        chart.tooltip.valueFormatter((y,index, p, d) => { 
            return 'TP='+y + ' Qual=' + fracFormat(d.point.quality) + ', Sens=' + percFormat(y / filteredTruth.length) + ' Prec='+percFormat(y / (d.point.x + y))  
        })
        
        chart.tooltip.headerFormatter(function(d) { 
            return 'False Positives = ' + d
        })
  
        
        console.log("rendering to " + id);
        d3.select('#' + id)
          .datum(points)  
          .call(chart);  
        
        let yScale = d3.scale.linear()
                             .domain(chart.lines.yScale().domain())
                             .range(chart.lines.yScale().range())
                            
        // Add right hand axis with percentage sensitivity
        var axis = nv.models.axis()
                            .scale(yScale)
                            .orient('right')
                            .tickPadding(6)
                            .tickValues([0,10,20,30,40,50,60,70,80,90,100].map(x => Math.round(filteredTruth.length*(x/100))))
                            .tickFormat(y => { console.log('tick y = ' + y); let val = percFormat(y / filteredTruth.length); if(y==filteredTruth.length) { return ' ' + val + ' Sensitivity';}; return val;})
                

        d3.select('#'+id+' .nv-wrap.nv-lineChart .nv-focus')
          .selectAll('.nv-y2')
          .data([points])      
          .enter()
          .append('g')
          .attr('class', 'nv-y2 nv-axis')
          .attr('transform', 'translate(' + (chart.xAxis.scale().range()[1]+10) + ',0)') 
          .call(axis);   
                
        
        window.chart = chart;
//        nv.utils.windowResize(function() { chart.update() });        
    }
}

Vue.component('roc-curve', {
    
    props: {
        model: { default: null },
        sliderwidth: { default: '180px' }
    },
    
    computed: {
        callerCount: function() {
           let count = Object.keys(this.model.cnv_calls).filter(c => c != 'truth').length
           return count
        }
    },
    
    watch: {
        targetRange: _.debounce(function() {
            this.renderPlot()
        },1000),
        
        callTargetRange: _.debounce(function() {
            this.renderPlot()
        },1000),
        
         sizeRange: _.debounce(function() {
            this.renderPlot()
        },1000),
        
        callSizeRange: _.debounce(function() {
            this.renderPlot()
        },1000)  
    },
    
    methods: {
       
        sizeValue: function(n) {
            return humanSize(Math.pow(10,n))    
        },
        
        targetStops: function(i) {
            let stops = [1,2,3,5,10,20,50,-1];    
            let stop = stops[i];
            if(stop < 0)
                return "Infinity";
            else
                return stop
        },
        
        formatSizeTooltip: function(x) {
            return humanSize(Math.pow(10,x[0])) + " - " + humanSize(Math.pow(10,x[1]))
        },
        
        formatTargetsTooltip: function(x) {
            return this.targetStops(x[0]) + " - " + this.targetStops(x[1])
        },
        
        renderPlot: function() {
            let plot = new CNVROCCurve({
                cnvs: this.cnv_calls,
                sizeRange: this.sizeRange,
                targetRange: [this.targetStops(this.targetRange[0]), this.targetStops(this.targetRange[1])],
                callSizeRange: this.callSizeRange,
                callTargetRange: [this.targetStops(this.callTargetRange[0]), this.targetStops(this.callTargetRange[1])],
            });
            plot.render('cnv_roc_curve'); 
        },
        
        showROCCurve: function() {
            console.log("ROC Curve, slider width = " + this.sliderwidth);
            
            $('#cnv_roc_curve_container').html(
              '<svg style="display: inline;" id=cnv_roc_curve></svg>' 
            );
            this.renderPlot()
        }
    },
    
    data: function() { return this.model },
    
    template: `
       <div>
        <h2>ROC Style Curves for {{callerCount}} CNV Callers</h2>
        <p>ROC curves show the number of true positives vs false positives found
           as quality score (or confidence measure) decreases.</p>
                       
        <div id=cnv_roc_curve_container>
            <svg id=cnv_roc_curve>
            </svg>
        </div>
          
            <div class=filterSettings>
            
            <n3-container fluid>
                 <n3-row>
                      <n3-column :col="12" class="header">
                          <h3>Filter Settings</h3>
                      </n3-column>
                 </n3-row>
                 
                 <n3-row>
                      <n3-column :col="6" class="context">
                          <h4>True Positives</h4>
                     </n3-column>
                      <n3-column :col="6" class="context">
                          <h4>CNV Calls</h4>
                     </n3-column>
                 </n3-row> 
                 <n3-row>
                      <n3-column :col="6" class="context">
                          <div id=newsizerange>
                              <div id=size_slider_label2 class='sliderLabel'>
                                 Size Range: {{sizeValue(sizeRange[0])}} - {{sizeValue(sizeRange[1])}} 
                              </div>
                              <n3-slider v-model="sizeRange" tooltip="always" :formatter='formatSizeTooltip' 
                                         :range="true" :min=0 :max=7 :width='sliderwidth'></n3-slider>
                          </div> 
                          <div id=newtargetsrange>
                              <div id=target_slider_label2 class='sliderLabel'>
                                  Target Regions: {{targetStops(targetRange[0])}}  - {{targetStops(targetRange[1])}}
                              </div>
                              <n3-slider v-model="targetRange" tooltip="always" :formatter='formatTargetsTooltip' 
                                         :range="true" :min=0 :max=7 :width='sliderwidth'></n3-slider>
                          </div>
                      </n3-column>
                      <n3-column :col="6">
                          <div id=newsizerange>
                              <div id=size_slider_label2 class='sliderLabel'>
                                 Size Range: {{sizeValue(callSizeRange[0])}} - {{sizeValue(callSizeRange[1])}} 
                              </div>
                              <n3-slider v-model="callSizeRange" tooltip="always" :formatter='formatSizeTooltip' 
                                        :range="true" :min=0 :max=7 :width='sliderwidth'></n3-slider>
                          </div> 
                          <div id=newtargetsrange>
                              <div id=target_slider_label2 class='sliderLabel'>
                                 Target Regions: {{targetStops(callTargetRange[0])}}  - {{targetStops(callTargetRange[1])}}
                              </div>
                              <n3-slider v-model="callTargetRange" tooltip="always" :formatter='formatTargetsTooltip' 
                                    :range="true" :min=0 :max=7 :width='sliderwidth'></n3-slider>
                          </div>                          
                      </n3-column>                
                  </n3-row>            
            </n3-container>
        </div>
      </div>
    `
})

