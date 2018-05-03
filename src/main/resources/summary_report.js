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
 * TODO: these should really get uploaded by the code since the 
 * reference could be non-hg19 compatible.
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
            cnv.size = cnv.end - cnv.start;
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
            window.model.cnv_calls = results;
            vue.cnv_calls = results
            window.rare_calls = {};
            
            // A number of different outputs are based on CNVs absent from the population,
            // so 
            Object.keys(cnv_calls).filter(caller => caller != 'truth').forEach((caller) => {
                window.rare_calls[caller] = cnv_calls[caller].filter(cnv => cnv.spanningFreq < MAX_RARE_CNV_FREQ)
            });
            
            $('.loading').remove()  
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
       if(!window.cnv_calls) {
           $('#' + activePanelId).prepend(`<div class='loading'><div class=loadingmsg>Loading, please wait!</div></div>`)
       }
       loadCnvs(fn, runs.map((r) => r)) // hack to clone runs
    }
    else {
        fn();
    }
}

function loadAndShowTab(id, callback) {
    loadAndCall(callback)
    window.location.hash = id
}

var activePanelId = null;

function activatePanel(panelId) {
    
   window.activePanelId = panelId;
       
   if(panelId == "qscorecalibration") {
       loadAndShowTab(panelId,showQscores);
   }
   else 
   if(panelId == "simrocs") {
       loadAndShowTab(panelId,vue.$refs.rocCurve.showROCCurve);
   }
   else 
   if(panelId == "sample_counts") {
       loadAndShowTab(panelId,vue.$refs.sampleCounts.showSampleCounts);
   }
   else 
   if(panelId == "genome_dist") {
       loadAndShowTab(panelId,function() { vue.$refs.genomeDistribution.showCNVGenomeDistribution() });
   } 
   else
   if(panelId == "sizebreakdown") {
       loadAndShowTab(panelId,showSizeBreakdown);
   }
   else
   if(panelId.match(/runcalls[0-9]*/)) {
       let runIndex = panelId.match(/runcalls([0-9])*/)[1];
           
       let runId = runs[parseInt(runIndex,10)];
           
       console.log(`Show calls ${runId}`);
           
       loadCnvReport(runId, () => showCNVReport(runIndex, runId));
   }
   else {
       console.log('Warning: unknown panel ' + panelId + ' activated')
   }
}

window.model = {
   analysisconfig: analysisConfig,
   cnv_calls: {},
   targetRange: [0,7],
   sizeRange: [0,7],
   callSizeRange: [0,7],
   callTargetRange: [0,7],
   maxFreq : MAX_RARE_CNV_FREQ,
   
   /**
    * Lists of excluded samples, keyed by caller id
    */
   excludedSamples: {},
   
   targetStops: function(i) {
        let stops = [1,2,3,5,10,20,50,-1];    
        let stop = stops[i];
        if(stop < 0)
            return "Infinity";
        else
            return stop
   },
 
   filterCNVs: function(simulationRegionsOnly) {
        let rawCnvs = this.cnv_calls;
        let callSizeMin = Math.pow(10, this.callSizeRange[0]);
        let callSizeMax = Math.pow(10, this.callSizeRange[1]);
        let callTargetsMin = this.targetStops(this.callTargetRange[0]);
        let callTargetsMax = this.targetStops(this.callTargetRange[1]);
        
        const unfilteredCount = Object.values(rawCnvs).reduce((n,caller) => n+caller.length, 0);
        console.log(`There are ${unfilteredCount} raw cnv calls`);
        
        console.log("Filtering by spanningFreq < " + this.maxFreq + " size range = " + this.sizeRange);
        
        this.filteredCnvs = {};
        
        let filters = [
            (cnv) => (cnv.targets >= callTargetsMin) && (cnv.targets<=callTargetsMax),
            (cnv) => (cnv.end - cnv.start > callSizeMin) && (cnv.end - cnv.start < callSizeMax),
        ]
        
        if(!Object.values(this.excludedSamples).every(sampleList => sampleList.length == 0))  {
            console.log('some samples are filtered')
            filters.push((cnv,caller) => !this.excludedSamples[caller] || this.excludedSamples[caller].indexOf(cnv.sample)<0)
        }
        
        if(simulationRegionsOnly) {
            let simRegionFilter = (cnv) => (cnv.chr != 'chrX' && cnv.chr != 'X')
            if(simulationType == 'replace') {
                simRegionFilter = (cnv) => (cnv.chr == 'chrX' || cnv.chr == 'X')
            }
            filters.push(simRegionFilter)
        }
        
        Object.keys(rawCnvs).filter(caller => caller != 'truth').forEach((caller) =>
            this.filteredCnvs[caller] = 
                rawCnvs[caller].filter(cnv => filters.every(fn => fn(cnv,caller)))
                               .sort((cnv1,cnv2) => cnv2.quality - cnv1.quality)
        );        

        return this.filteredCnvs
   }
}

var vue = null

main_template = `
<div>
   <div id=simrocs>
        <roc-curve :model='model' ref='rocCurve'/>
    </div><!-- simrocs -->
                    
    <div id=sample_counts>
        <sample-counts :model='model' ref='sampleCounts'/>
    </div> <!-- sample_counts -->
                    
    <div id=genome_dist>
        <genome-distribution :model='model' ref='genomeDistribution'/>
    </div>
                  
    <div id=sizebreakdown>
        <sens-by-size :model='model'/>
    </div> <!-- sizebreakdown -->
                    
    <div id=qscorecalibration>
        <qual-score-calibration :model='model'/>
    </div> <!-- qscore calibration --> 
</div>
`
    
function createTabs() {
    console.log("initializing tabs")
    window.tabs = $('#innerLayout').tabs({
            activate: function(event,ui) {
                                    
                console.log("Activate tab");
                $(events).trigger("activate",ui);
                                    
                if(!ui.newPanel[0]) {
                    console.log("No new panel");
                    return;
                }
                                    
                $('#maincenter').css('padding','0px');
    
                var panelId = ui.newPanel[0].id;
            }
    });
    
    // Ugly hack: inline width setting causes incorrect width - remove it
    $('.ui-layout-center').css('width', '')
}

require(['jquery.layout.min', 'N3Components.min','vue','lodash.min', 'd3', 'nv.d3','roccurve','ximmer_qc'], 
        function(jquery_layout, N3Components, Vue, lodash, d3js, nvd3, roccurve, ximmer_qc) {
    
    window.Vue = Vue;
    
    window.showSizeBreakdown = ximmer_qc.showSizeBreakdown
    window.showQscores = ximmer_qc.showQscores
    
    $(document).ready(function() {
        
       $('body').layout({ applyDefaultStyles: true });
       layout = $('#innerLayout').layout({ applyDefaultStyles: true });
       layout.sizePane("north",50); 
        
       console.log('Summary report init');
       
       console.log('Creating main vue: ' + Vue.hey)
       vue = new Vue({
           el: '#cnvsummarytabs',
           
           template: main_template,
           
           mounted: function() {
               console.log('main template rendered')
               
                // Calling tabs() directly doesn't seem to work
                setTimeout(createTabs, 0);
           },
           
           data: {
               model: model
           }
       }) 
    
       if((window.location.hash != '') && (window.location.hash != '#')) {
           let panelId = window.location.hash.replace(/^#/,'')
           console.log('Loading panel ' + panelId + ' from hash')
           activatePanel(panelId)
       }
       
       $(events).on("activate", (event,ui) => {
           
           console.log("Activated in summary");
           var panelId = ui.newPanel[0].id;
           activatePanel(panelId)
       }); 
    })
 })
