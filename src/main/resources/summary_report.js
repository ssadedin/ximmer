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

function loadCnvs(callback) {
    
    console.log("loadCnvs");
    
     var oldScriptElement = document.getElementById('cnvs_load_script');
     if(oldScriptElement)
         oldScriptElement.parentNode.removeChild(oldScriptElement);

     const script = document.createElement("script");
     script.id = 'cnvs_load_script';
     script.src = '1/analysis-base/report/cnv_calls.js'
     script.async = true;
     script.onload = callback;
     
     console.log("Loading cnv calls from " + script.src);
     document.body.appendChild(script); 
}

function loadAndShowQscores() {
    console.log("loadAndShowQscores");
    
    if(typeof(window.cnv_calls) == 'undefined') {
        loadCnvs(showQscores)
    }
    else {
        showQscores();
    }
}

$(document).ready(function() {
   console.log('Summary report init');
   $(events).on("activate", (event,ui) => {
       console.log("Activated in summary");
       var panelId = ui.newPanel[0].id;
       if(panelId == "qscorecalibration") {
           loadAndShowQscores()
       }
   }); 
})