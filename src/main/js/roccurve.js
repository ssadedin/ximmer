
// ----------------------------- ROC Curve Functions  ----------------------------------

define(["vue","nv.d3","N3Components.min","genomicranges"], function(Vue, nvd3, N3Components,genomicranges) {
    
    console.log("Register N3:")
    console.log(N3Components)
    Vue.use(N3Components.default, 'en');

    class CNVROCCurve {
        
        constructor(props) {
            /**
             * cnvs - a map keyed on caller id with values being the calls
             */
            this.rawCnvs = props.cnvs;
            if(!this.rawCnvs)
                return
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
        
        calculateCombinations(cnvs) {
            
            let union = []
            let intersect = []
            
            let cfgLabels = window.model.analysisConfig.callerCfgs
            let intersectCallerWeights = cfgLabels.reduce((acc,c) => { 
                acc[c] = model.rocIntersectCallers.indexOf(c)>=0?1:0
                return acc;
             }, {})
                
            let unionCallerWeights = cfgLabels.reduce((acc,c) => {
                acc[c] = model.rocUnionCallers.indexOf(c)>=0?1:0
                return acc;
            }, {})            
            
            // Clear any previous flags, and set ranks and caller ids
            
            const callerIdMap = [
                [/^conifer_/, 'cfr_'],
                [/^exomedepth_/,'ed_'],
            ]
            
            let cnvLists = Object.values(cnvs)
            let callerIds = Object.keys(cnvs)
            
            const unionCallers = model.rocUnionCallers.map(callerId => {
                for(const mapping of callerIdMap) {
                    callerId = callerId.replace(mapping[0],mapping[1])
                }
                return callerId
            })
            
            const intersectCallers = model.rocIntersectCallers.map(callerId => {
                for(const mapping of callerIdMap) {
                    callerId = callerId.replace(mapping[0],mapping[1])
                }
                return callerId
            }) 
            
            cnvLists.forEach((callerCnvList, callerIndex) => {
                const callerId = callerIds[callerIndex]
                callerCnvList.forEach((cnv,rank) =>  { 
                    cnv.unioned = false; 
                    cnv.intersected=false 
                    cnv.rank = rank
                    cnv.caller = callerId
                })
            })
  
            // Create combined list of all callers,
            let allCnvs = [].concat.apply([], cnvLists)
            
            // Index the regions
            const cnvRegions = new Regions(allCnvs)
            
            // Sort by rank
            allCnvs.sort((cnv1, cnv2) => cnv1.rank - cnv2.rank)
            
            // Now iterate the whole list
            for(const cnv of allCnvs) {
                
                if(cnv.unioned && cnv.intersected)
                    continue
                    
                const overlaps = cnvRegions.getOverlaps(cnv).filter(otherCnv => cnv.sample == otherCnv.sample)
                
                if(!cnv.unioned && unionCallers.indexOf(cnv.caller)>=0) {
                    
                    // At least one of the union callers must have called it 
                    if(unionCallers.find(caller => overlaps.find(o => o.caller == caller))) {
                        overlaps.forEach(o => o.unioned = true)
                        // todo: expand range of what is pushed to min / max boundaries
                        union.push(Object.assign({},cnv))
                    }
                }
                
                if(!cnv.intersected && intersectCallers.indexOf(cnv.caller)>=0) {
                    if(intersectCallers.every(c => overlaps.find(o => o.caller == c))) {
                        overlaps.forEach(o => o.intersected = true)
                        intersect.push(Object.assign({},cnv))
                    }
                }
            }
            
            return {
                union: union,
                intersect: intersect
            }
        } 
        
        computeROCStats(cnvs,truthSet) {
            
            // the set of true positives that we have identified
            let tps = truthSet.map(function(cnv, i) {
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
                        // This case occurs when the tp set is filtered (eg: by size etc.)
                    }
                    else
                    if(!tp.detected) {
                        tp.detected = true;
                        ++tpCount;
                    }
                }
                else {
                    if(cnv.DGVFreq < MAX_RARE_CNV_FREQ)
                        ++fpCount;
                }
                cnv.tp = tpCount;
                cnv.fp = fpCount
            });
        }
        
        /**
         * Return a map of id to label
         */
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
            
            if(!this.rawCnvs)
                return
            
            this.filteredCnvs = model.filterCNVs(true)
            
            let rawCnvs = this.rawCnvs;
            let sizeMin = Math.pow(10, this.sizeRange[0]);
            let sizeMax = Math.pow(10, this.sizeRange[1]);
            let targetsMin = this.targetRange[0];
            let targetsMax = this.targetRange[1];
            
            // Note: heuristic here: if none of the true positives are on chrX then assume it was excluded
            // X chr might be excluded if a mix of males and females was used, since this would cause
            // an inflation of false positives 
            let includeChrX = (simulationType == 'replace') || 
                              !this.rawCnvs.truth.every(cnv => (cnv.chr != 'chrX') && (cnv.chr != 'X')) // no tp on chrX
            
            console.log("Include chrX? " + includeChrX)
            
            let filteredTruth = 
                this.rawCnvs.truth.filter((cnv) => (cnv.targets >= targetsMin) && (cnv.targets<=targetsMax) && 
                                                   (cnv.end - cnv.start > sizeMin) && (cnv.end - cnv.start < sizeMax) &&
                                                   (includeChrX || (cnv.chr != 'chrX'))
                                         )
            
            console.log(`There are ${filteredTruth.length} true cnv calls after filtering by size ${targetsMin}-${targetsMax} and simulation type chr (${simulationType})`);
            
            let cnvCount = Object.values(this.filteredCnvs).reduce((n,caller) => n+caller.length, 0);
            console.log(`There are ${cnvCount} cnv calls after filtering by spanningFreq<${this.maxFreq}`);
            
            window.cnvs = this.filteredCnvs;
                
            // Now iterate through each caller's CNVs and compute the number of true and false positives
            Object.values(this.filteredCnvs).forEach((cnvs) => this.computeROCStats(cnvs,filteredTruth));
            
            if(model.rocComboOptions.length>0) {
                console.log('Calculating combinations')
                let combos = 
                       this.calculateCombinations(this.filteredCnvs)
                       
                if(model.rocComboOptions.indexOf('union')>=0) {
                    this.filteredCnvs.union = combos.union
                    this.computeROCStats(this.filteredCnvs.union,filteredTruth)
                }
                if(model.rocComboOptions.indexOf('intersection')>=0) {
                    this.filteredCnvs.intersect = combos.intersect
                    this.computeROCStats(this.filteredCnvs.intersect,filteredTruth)
                }
                    
            }
            
            let callerLabels = this.computeCallerLabelMap()
           
            let points = [];
            Object.keys(this.filteredCnvs).forEach(caller => points.push({
                values: this.filteredCnvs[caller].map(cnv => { return { x: cnv.fp, y: cnv.tp, quality: cnv.quality }}),
                key: callerLabels.get(caller.replace(/\./,"_"))  || caller            
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
                    
            
            this.points = points;
            
            window.chart = chart;
    //        nv.utils.windowResize(function() { chart.update() });        
        }
    }
    
    console.log("Registering roc-curve component")
    
    Vue.component('roc-curve', {
        
        props: {
            model: { default: null },
            sliderwidth: { default: '180px' }
        },
        
        computed: {
            callerCount: function() {
               let count = Object.keys(this.model.cnv_calls).filter(c => c != 'truth').length
               return count
            },
            
            hiddenSamples: function() {
                let allHidden = {}
                for(var ex of Object.values(model.excludedSamples)) {
                    ex.forEach(s => allHidden[s]=1)
                }
                return Object.keys(allHidden)
            },
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
            },1000),
            
            rocComboOptions: _.debounce(function() {
                this.renderPlot()
            },1000),
            
            rocIntersectCallers: _.debounce(function() {
                this.renderPlot()
            },1000) ,
            
            rocUnionCallers: _.debounce(function() {
                this.renderPlot()
            },1000)  
        },
        
        methods: {
            
            downloadData: function() {
                var json = JSON.stringify({
                    truth: model.cnv_calls.truth, 
                    points: this.plot.points
                })
                if (json == null) return;
                if (!json.match(/^data:text\/json/i)) {
                    json = 'data:text/json;charset=utf-8,' + json;
                }
                let data = encodeURI(json);
                $('#downloadROCExport').attr('download', 'roccurve.json')
                                       .attr('href', json)
                $('#downloadROCExport')[0].click()
            },
           
            sizeValue: function(n) {
                return humanSize(Math.pow(10,n))    
            },
            
            formatSizeTooltip: function(x) {
                return this.sizeValue(x[0]) + " - " + this.sizeValue(x[1])
            },
            
            formatTargetsTooltip: function(x) {
                return this.targetStops(x[0]) + " - " + this.targetStops(x[1])
            
            },
            
            renderPlot: function() {
                let plot = new CNVROCCurve({
                    cnvs: this.cnv_calls,
                    sizeRange: this.sizeRange,
                    targetRange: [this.model.targetStops(this.targetRange[0]), this.model.targetStops(this.targetRange[1])],
                    callSizeRange: this.callSizeRange,
                    callTargetRange: [this.model.targetStops(this.callTargetRange[0]), this.model.targetStops(this.callTargetRange[1])],
                });
                plot.render('cnv_roc_curve'); 
                
                this.plot = plot;
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
            <span id='downloadLinkBox'>
                <a href='#' id=downloadROC v-on:click.prevent='downloadData()'>download</a>
                <a href='#' id='downloadROCExport'>dummy</a>
            </span>
           
            <h2>ROC Style Curves for {{callerCount}} CNV Callers</h2>
            
            <p v-if='hiddenSamples.length==0'>ROC curves show the number of true positives vs false positives found
               as quality score (or confidence measure) decreases.</p>
               
            <p v-else><b>Note:</b> these samples are hidden in one or more callers: {{hiddenSamples.join(',')}}</p>
                           
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
                          <n3-column :col="4" class="context">
                              <h4>True Positives</h4>
                         </n3-column>
                          <n3-column :col="4" class="context">
                              <h4>CNV Calls</h4>
                         </n3-column>
                         <n3-column :col="4" class="context">
                              <h4>Combinations</h4>
                         </n3-column> 
                     </n3-row> 
                     <n3-row>
                          <n3-column :col="4" class="context">
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
                          <n3-column :col="4">
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
                          <n3-column :col="4">
                              <n3-checkbox-group v-model="rocComboOptions">
                                <div>
                                    <!-- union -->
                                    <n3-checkbox class='combx_checkbox' label="union">Union</n3-checkbox>
                                    <n3-checkbox-group v-model="rocUnionCallers" v-if='rocComboOptions.indexOf("union")>=0'>
                                        <div class='subCheckList'>
                                            <n3-checkbox v-for='(caller,index) in analysisConfig.callerLabels' class='combx_checkbox' :label="analysisConfig.callerCfgs[index]">{{caller}}</n3-checkbox>
                                        </div>
                                    </n3-checkbox-group>  
                                    
                                    <!-- intersect -->
                                    <n3-checkbox class='combx_checkbox' label="intersection">Intersection</n3-checkbox>
                                    <n3-checkbox-group v-model="rocIntersectCallers" v-if='rocComboOptions.indexOf("intersection")>=0'>
                                        <div class='subCheckList'>
                                            <n3-checkbox v-for='(caller,index) in analysisConfig.callerLabels' class='combx_checkbox' :label="analysisConfig.callerCfgs[index]">{{caller}}</n3-checkbox>
                                        </div>
                                  </n3-checkbox-group> 
                                </div>
                              </n3-checkbox-group>
                         </n3-column>  
                      </n3-row>            
                </n3-container>
                
            </div>
        
          </div>
        `
    })
    
    console.log("Registered ROC Curve")

    return {
        CNVROCCurve : "CNVROCCurve"
    }
})
    
