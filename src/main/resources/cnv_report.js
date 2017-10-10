// vim: ts=4:sw=4:cindent:expandtab
/*!
 * jQuery create HTML plugin
 * Original author: @efunneko
 * Licensed under the MIT license
 */


// Simple jQuery plugin to make it a bit cleaner to add
// HTML using jQuery. Instead of doing:
//
// $("body").$("<div/>").append(
//   $("<span/>", {id: 'mySpan'}).append(
//     $("<b/>").html("My bold statement!")));
//
// You can do:
//
// $("body").$div().$span({id: 'mySpan'}).b("My bold statement!");

(function ( $, window, document, undefined ) {

    // Create the defaults once
    var pluginName = 'createHtml',
        defaults = {
        },
        elements = [
            "a", "abbr", "address", "area", "article", "aside", "audio", "b", "base",
            "bdi", "bdo", "blockquote", "body", "br", "button", "canvas", "caption",
            "cite", "code", "col", "colgroup", "command", "data", "datalist", "dd",
            "del", "details", "dfn", "div", "dl", "dt", "em", "embed", "fieldset",
            "figcaption", "figure", "footer", "form", "h1", "h2", "h3", "h4", "h5",
            "h6", "head", "header", "hgroup", "hr", "html", "i", "iframe", "img", "input",
            "ins", "kbd", "keygen", "label", "legend", "li", "link", "main", "map", "mark", "math",
            "menu", "meta", "meter", "nav", "noscript", "object", "ol", "optgroup", "option",
            "output", "p", "param", "pre", "progress", "q", "rp", "rt", "ruby", "s",
            "samp", "script", "section", "select", "small", "source", "span", "strong",
            "style", "sub", "summary", "sup", "svg", "table", "tbody", "td", "textarea",
            "tfoot", "th", "thead", "time", "title", "tr", "track", "u", "ul", "var",
            "video", "wbr"
        ],
        methods = {
            configure: function(options) {
                if (options.installParentFunctions) {
                    // Add $<element>_ functions that will return the parent object 
                    // rather than itself
                    $.each(elements, function(i, elName) {
                        $.fn["$" + elName + "_"] = function(content, attrs) {
                            if (typeof(content) == 'object' && typeof(attrs) == 'undefined') {
                                attrs = content;
                                content = undefined;
                            }
                            var el = $("<" + elName + ">", attrs);
                            if (content && content != "") {
                                el.html(content);
                            }
                            var parent = $(this);
                            parent.append(el);
                            return parent;
                        }
                    });
                }
            }
        };

    $[pluginName] = function(method, options) {
        if (methods[method]) {
            methods[method](options);
        }
    };

    // Add all element functions
    $.each(elements, function(i, elName) {
        $.fn["$" + elName] = function(content, attrs) {
            if (typeof(content) == 'object' && typeof(attrs) == 'undefined') {
                attrs = content;
                content = undefined;
            }
            var el = $("<" + elName + ">", attrs);
            if (content && content != "") {
                el.html(content);
            }
            $(this).append(el);
            return el;
        }
    });
    
})( jQuery, window, document );

function partial(fn) {
    var args = Array.prototype.slice.call(arguments);
    args.shift();
    return function() {
        var new_args = Array.prototype.slice.call(arguments);
        args = args.concat(new_args);
        return fn.apply(window, args);
    };
}

/*
 * jqModal - Minimalist Modaling with jQuery
 *
 * Copyright (c) 2007-2015 Brice Burgess @IceburgBrice
 * Dual licensed under the MIT and GPL licenses:
 *   http://www.opensource.org/licenses/mit-license.php
 *   http://www.gnu.org/licenses/gpl.html
 * 
 * $Version: 1.3.0 (2015.04.15 +r24)
 * Requires: jQuery 1.2.3+
 */

!function(n){n.fn.jqm=function(t){return this.each(function(){var o=n(this),e=o.data("jqm")||n.extend({ID:m++},n.jqm.params),i=n.extend(e,t);o.data("jqm",i).addClass("jqm-init")[0]._jqmID=i.ID,i.trigger&&o.jqmAddTrigger(i.trigger)})},n.fn.jqmAddTrigger=function(o){return this.each(function(){s(n(this),"jqmShow",o)||t("jqmAddTrigger must be called on initialized modals")})},n.fn.jqmAddClose=function(o){return this.each(function(){s(n(this),"jqmHide",o)||t("jqmAddClose must be called on initialized modals")})},n.fn.jqmShow=function(t){return this.each(function(){!this._jqmShown&&o(n(this),t)})},n.fn.jqmHide=function(t){return this.each(function(){this._jqmShown&&e(n(this),t)})};var t=function(n){window.console&&window.console.error&&window.console.error(n)},o=function(t,o){var e=t.data("jqm"),o=o||window.event,i=parseInt(t.css("z-index")),i=i>0?i:3e3,a=n("<div></div>").addClass(e.overlayClass).css({height:"100%",width:"100%",position:"fixed",left:0,top:0,"z-index":i-1,opacity:e.overlay/100}),s={w:t,c:e,o:a,t:o};if(t.css("z-index",i),e.ajax){var d=e.target||t,c=e.ajax;d="string"==typeof d?n(d,t):n(d),"@"==c.substr(0,1)&&(c=n(o).attr(c.substring(1))),d.load(c,function(){e.onLoad&&e.onLoad.call(this,s)}),e.ajaxText&&d.html(e.ajaxText),r(s)}else r(s)},e=function(n,t){var o=n.data("jqm"),t=t||window.event,e={w:n,c:o,o:n.data("jqmv"),t:t};d(e)},i=function(t){return t.c.overlay>0&&t.o.prependTo("body"),t.w.show(),n.jqm.focusFunc(t.w,!0),!0},a=function(n){return n.w.hide()&&n.o&&n.o.remove(),!0},s=function(t,o,e){var i=t.data("jqm");return t.data("jqm")?n(e).each(function(){this[o]=this[o]||[],n.inArray(i.ID,this[o])<0&&(this[o].push(i.ID),n(this).click(function(){return t[o](this),!1}))}):!1},r=function(t){var o=t.w,e=t.o,i=t.c;i.onShow(t)!==!1&&(o[0]._jqmShown=!0,i.modal?(!l[0]&&c("bind"),l.push(o[0])):o.jqmAddClose(e),i.closeClass&&o.jqmAddClose(n("."+i.closeClass,o)),i.toTop&&e&&o.before('<span id="jqmP'+i.ID+'"></span>').insertAfter(e),o.data("jqmv",e),o.unbind("keydown",n.jqm.closeOnEscFunc),i.closeOnEsc&&o.attr("tabindex",0).bind("keydown",n.jqm.closeOnEscFunc).focus())},d=function(t){var o=t.w,e=t.o,i=t.c;i.onHide(t)!==!1&&(o[0]._jqmShown=!1,i.modal&&(l.pop(),!l[0]&&c("unbind")),i.toTop&&e&&n("#jqmP"+i.ID).after(o).remove())},c=function(t){n(document)[t]("keypress keydown mousedown",u)},u=function(t){var o=n(t.target).data("jqm")||n(t.target).parents(".jqm-init:first").data("jqm"),e=l[l.length-1];return o&&o.ID==e._jqmID?!0:n.jqm.focusFunc(e,t)},m=0,l=[];n.jqm={params:{overlay:50,overlayClass:"jqmOverlay",closeClass:"jqmClose",closeOnEsc:!1,trigger:".jqModal",ajax:!1,target:!1,ajaxText:"",modal:!1,toTop:!1,onShow:i,onHide:a,onLoad:!1},focusFunc:function(t,o){return o&&n(":input:visible:first",t).focus(),!1},closeOnEscFunc:function(t){return 27==t.keyCode?(n(this).jqmHide(),!1):void 0}}}(jQuery);

var initialized = false;

var cnvTable = null;
var filteredTable = null;
var cnvData = [];
var filteredData = [];
var columns = [];

// Array of strings to 'eval' as filters
var filters = [];
var layout = null;
var cnvLayout = null;
var cnvIndex = -1;


var geneList = {};

var allTags = [];
var userAnnotations = { 
    tags : {} 
};

var savedSettings = {};

if(localStorage.savedCnvSettings) {
   savedSettings = JSON.parse(localStorage.savedCnvSettings);
}

if(!localStorage.cnvSettings) {
    var defaultSettings = { }
    defaultSettings[location.href] = { filters: [] }
    localStorage.cnvSettings = JSON.stringify(defaultSettings);
}

var SAMPLE_ID_COLUMN=7;

var idMaskRegExp = (typeof(window.idMask)!='undefined') ? new RegExp(idMask) : null;

function createCnvRow(row, data, dataIndex) {
    
    var tds = row.getElementsByTagName('td');

    var index = parseInt(row.getElementsByTagName('a')[0].innerHTML,10)-1;
    
    var cnvTags = userAnnotations[index] && userAnnotations[index].tags;
    addRowTags(row, cnvTags);
    if(idMaskRegExp) {
        var match = data[SAMPLE_ID_COLUMN].match(idMaskRegExp);
        if(match)
            tds[SAMPLE_ID_COLUMN].innerHTML = match[1];
    }

    var cnv = cnvs[index];
    
    // highlight row if overlapping truncating variant
    if(cnv.truncated) {
        $(row).addClass('truncated')
    }
    
    for(key in geneList) {
        if(cnv.genes.indexOf(key)>=0) {
            $(row).addClass('genelist')
            $(row).addClass('genelist'+geneList[key])
        }
    }
}

function addRowTags(row, cnvTags) {
    var td0 = row.getElementsByTagName('td')[0];
    
    if(typeof(td0.oldHTML) == 'undefined')
        td0.oldHTML = td0.innerHTML;
    
    var tagDivs = [];
    for(tag in cnvTags) {
        tagDivs.push('<span class="rowTagDiv tag'+allTags.indexOf(tag)+'">' + tag + '</span>')
    }
    td0.innerHTML = tagDivs.join(' ') + td0.oldHTML;
}

//////////////////////////////////////// VueJS Components ///////////////////////////////////////

/*
Vue.component('CNVTags', {
    
    template: `
        
    `
})
*/



//////////////////////////////////////// Main Entry Point ///////////////////////////////////////

$(document).ready(function() {

    if(initialized)
        return;
    
    if(typeof(window.cnvs) == 'undefined') {
        console.log('CNVs not defined for page: do not set up CNV report');
        return;
    }
    
    console.log("init");

    columns = []
    $('#cnvTable thead th').each(function() {columns.push(this.innerHTML);})
    columns = columns.map(function(x) { return {title: x}; })

    cnvTable = $('#cnvTable').DataTable({ 
        "iDisplayLength": 25,
        "columnDefs": [
            { "width": "2em", "targets": 0 },
            { "width": "25%", "targets": 0 },
            { "width": "6em", "targets": 0 },
            { "width": "3em", "targets": 0 }
        ]
    });

    cnvData = cnvTable.rows().data();
    
    $('body').layout({ applyDefaultStyles: true });

    cnvLayout = layout = $('#innerLayout').layout({ applyDefaultStyles: true });
    layout.sizePane("north",50);
    
    $(".cnvimg").click(function() {
        window.open(this.src, "cnvimg", "width=1600,height=1024");
    });

    $('#filterHVR').change(function() {
        showHVR(this.checked);
    });

    add_display_events();

    renderFilters();

    with($('#filterOuter')) {
        $button({id:'addFilter'}).$span('Add a filter');
        $button({id:'runFilters'}).$span('Run');
        with($button({id:'clearFilters'})) {
            $span('Clear');
            click(function() { filters=[]; renderFilters(); filterTable(); $('#filterHVR').removeAttr('checked'); });
        }
        with($button({id:'clearTags'})) {
            $span('Clear Tags');
            click(function() { for(i in userAnnotations) { userAnnotations[i].tags = {}; }; filterTable(); });
        } 
        $button('Share').click(function() {
            $('#message').html('');
            with($('#message').$div()) {
                $h4('Copy the text below to share these settings with somebody else:');
                $button('Copy').click(function() {
                    var range = document.createRange();  
                    range.selectNode($('#copyDiv')[0]);  
                    window.getSelection().addRange(range); 
                    if(!document.execCommand('copy')) {
                        alert('Sorry, could not copy the data to your clipboard');
                    }
                });
                $br();
                $div({id:'copyDiv'}).$span(JSON.stringify(JSON.parse(localStorage.cnvSettings)[location.href]));
            }
            $('#message').jqm();
            $('#message').jqmShow();
        });
        $button('Import').click(function() {
            var newSettings = prompt("Paste filters from another session below:");
            if(newSettings) {
                loadSettings(JSON.parse(newSettings));
                saveSettings();
            }
        });
        $button('Gene List').click(function() {
            
            var formatted = [];
            var catToGenes = {
            }
            for(gene in geneList) {
                var cat = geneList[gene];
                if(typeof(catToGenes[cat]) == 'undefined')
                    catToGenes[cat] = [];
                catToGenes[cat].push(gene);
            }

            for(var cat in catToGenes) {
                formatted.push(cat + ':' + catToGenes[cat].join(","));
            }

            $('#genelist')[0].value = formatted.join('\n');
            
            $('#dialog').jqm();
            $('#dialog').jqmShow();
        });

        $button('Save As').click(function() {
            var saveName = prompt('Please enter a name to save your filters and tags with: ','');
            savedSettings[saveName] = JSON.parse(localStorage.cnvSettings)[location.href];
            localStorage.savedCnvSettings = JSON.stringify(savedSettings);
            renderSavedSettings();
        });
        $span({id:'savedSettingsOuter'});
        renderSavedSettings();
        
        $button('Help').click(function() {
            $('#help').jqm();
            $('#help').jqmShow();
        });
    }

    $('#geneOk').click(function() {
        $('#dialog').jqmHide();
        geneList = {}
        var values = $('#genelist')[0].value.split(/[ \n\r]/).map(function(e) {
            var geneSplit = e.split(':');
            if(geneSplit.length > 1) {
                var cat = parseInt(geneSplit[0],10);
                var genes = geneSplit[1].split(',');
                genes.forEach(function(gene) {
                    geneList[gene] = cat;
                });
            }
            else {
                geneList[geneSplit[0]] = 1;
            }
            saveSettings();
            updateGeneList();
            filterTable();
        });

    });

    $('#geneCancel').click(function() {
        $('#dialog').jqmHide();
    });

    $('#addFilter').click(function() {
        filters.push({ expr: '', id: filters.length});
        renderFilters();
    });
    $('#runFilters').click(filterTable);

    var attrs = $.map(cnvs[0], function(v,k) { return k=='chr' ? 'chr' : k });
    $.each(annotations[0],function(k,v) { attrs.push(k); }); 

    $('#addFilter').one("click",function() {
        with($('#filterHelp')) {
            $span("Filter attributes: " + attrs.join(","));
        }
        layout.sizePane("north",layout.panes.north.outerHeight() + 10);
        $('#filterHelp').slideDown();
    });

    if(localStorage.cnvSettings) {
        console.log("Load settings for: " + location.href);
        var oldSettings = JSON.parse(localStorage.cnvSettings)[location.href];
        if(oldSettings) {
            loadSettings(oldSettings);
        }
    }
    
    var cnvIndex = 1;
    $('#cnvTable').on('order.dt',  function() { setTimeout(add_display_events,0);});
    $('#cnvTable').on('page.dt',  function() { setTimeout(add_display_events,0);});
    $('#cnvTable').on('search.dt',  function() { setTimeout(add_display_events,0);});
    
    $(document.body).keydown(function(e) {
       if(e.keyCode == 74) { // j
           console.log("Key pressed: " + e.keyCode);
           if(highlightTr) {
               // Find index of the highlighted row in the table
               var trs = $('#cnvTable tr');
               for(var i=0; i<trs.length; ++i) {
                   if(trs[i] == highlightTr) {
                       // Find the next row
                       if(i<trs.length-1) {
                           console.log("go to next CNV: " + trs[i+1]);
                           showCNVDetailsForRow(trs[i+1]);
                       }
                       break;
                   }
               }
           }
       }
       else
       if(e.keyCode == 75) { // k
           
           if(highlightTr) {
               // Find index of the highlighted row in the table
               var trs = $('#cnvTable tr');
               for(var i=trs.length-1; i>=0; --i) {
                   if(trs[i] == highlightTr) {
                       // Find the next row
                       if(i>0) {
                           console.log("go to prev CNV: " + trs[i-1]);
                           showCNVDetailsForRow(trs[i-1]);
                       }
                       break;
                   }
               }
           } 

          console.log("Key pressed: " + e.keyCode); 
       }
    });

    initialized = true;
//    $('#tableHolder')[0].style.display='block';
});

//////////////////////////////////////// Subroutines ///////////////////////////////////////

function renderSavedSettings() {
    $('#savedSettingsOuter').html('');
    if(Object.keys(savedSettings).length==0)
        return;
    with($('#savedSettingsOuter')) {
        with($select()) {
          $option('Select Saved ...');
          for(n in savedSettings) {
            var o = $option({value:n});
            o.$span(n);
          }
          change(function() {
              var settingsName = this.value;
              var settings = savedSettings[settingsName];
              loadSettings(settings);
              saveSettings();
              $('#deleteSaved').html('');
              $('#deleteSaved').$button('Delete ' + settingsName).click(function() {
                if(confirm('Delete saved settings ' + settingsName + '?')) {
                    delete savedSettings[settingsName];
                    localStorage.savedCnvSettings = JSON.stringify(savedSettings);
                    renderSavedSettings();
                }
              });
          });
        }
        $span({id:'deleteSaved'});
    }
}

function loadSettings(oldSettings) {
    console.log("Restoring settings ...");
    filters = oldSettings.filters;
    geneList = oldSettings.geneList || {};
    
    if(oldSettings.userAnnotations) {
        userAnnotations = oldSettings.userAnnotations;
        for(var cnvIndex in userAnnotations) {
            for(var tag in userAnnotations[cnvIndex].tags) {
                if(allTags.indexOf(tag)<0)
                    allTags.push(tag);
            }
        }
    }
    console.log("User annotations are: " + userAnnotations);
    renderFilters();
    updateGeneList();
    filterTable();
}

function updateGeneList() {
    for(var i=0; i<cnvs.length; ++i) {
        cnvs[i].category = Math.max.apply(Math, cnvs[i].genes.map(function(gene) { return geneList[gene] ? geneList[gene] : 0; } ) )
    }
}

function toggleAmplicons() {
    var img = $('#coveragePlot img')[0];
    if(img.src.match(/\.ac\.png$/)) { // amplicon mode enabled
        img.src = img.src.replace(/\.ac\.png$/, '.png');
        $('#ampliconButton').html('Amplicons');
    }
    else {
        img.src = img.src.replace(/\.png$/, '.ac.png');
        $('#ampliconButton').html('Amplicons off');
    }
}

function renderFilters() {
    $('#filters').html('');
    with($('#filters').$span()) {
        for(var i=0; i<filters.length; ++i) {
            var f = filters[i];
            if(typeof(f.id) == 'undefined')
                continue
            $input({type:'text', id: f.id, value: f.expr}).keyup(function(e) {
                    if(e.keyCode == 13)
                        filterTable();
            }).focus();
        }
    }
}

function saveSettings() {
    for(var i=0; i<filters.length; ++i) {
        if(filters[i].id != null)
            filters[i].expr = document.getElementById(filters[i].id).value;
    }

    var cnvSettings = JSON.parse(localStorage.cnvSettings);
    if(!cnvSettings[location.href])
        cnvSettings[location.href] = {};
    
    cnvSettings[location.href].filters = filters;
    cnvSettings[location.href].geneList = geneList;
    cnvSettings[location.href].userAnnotations = userAnnotations;
    localStorage.cnvSettings = JSON.stringify(cnvSettings);
}


var dataMap = {
    'TRUE' : true,
    'FALSE' : false
}

function filterTable() {

    saveSettings();
    
    // It's irritating if the text in search box disappears when adjusting other filters
    var oldSearchText = $('#cnvTable_filter input')[0].value;

    // Update the filters from the fields
    // Filter out rows
    var newTable = cnvData.filter(function(value,index) { 
        var cnv = cnvs[index];
        cnv.variants = cnvVariants[index]
        cnv.index = index+1;
        window.dataval = value;
        var data = $.extend({}, cnv, annotations[index]);
        
        var cnvTags = userAnnotations[index] && userAnnotations[index].tags;
        data.tags = cnvTags ? cnvTags : {};
        
        data.size = (data.end - data.start);

        if(data.spanning == null)
            data.spanning = [];

        for(x in data) {
          if(typeof(dataMap[data[x]]) != 'undefined')
            data[x] = dataMap[data[x]];
        }

        var indexp = index+1;
        $('#cnv_'+indexp+'_detail').show();
        $('#cnv_'+indexp+'_img').show();

        with(data) {
            for(var i=0; i<filters.length; ++i) {
              var result;
              if(filters[i].expr == '')
                result = true;
              else {
                try { result = eval(filters[i].expr); }
                catch(e) { result = false; }
              }

              if(!result) {
                $('#cnv_'+indexp+'_detail').hide();
                $('#cnv_'+indexp+'_img').hide();
                return false;
              }
            }
        }
        return true
    });

    $('#tableHolder').html('<table id=cnvTable class=stripe></table>')
    filteredCnvTable = $('#cnvTable').DataTable({ data: newTable, 
                               columns: columns, 
                               iDisplayLength: 25,
                               destroy: true,
                               createdRow: createCnvRow
                               });
    filteredTable = cnvTable.rows().data();
    
    if(oldSearchText) {
        filteredCnvTable.search(oldSearchText).draw()
        $('#cnvTable_filter input')[0].value = oldSearchText
    }
    
    add_display_events();

    $('#cnvTable').on('order.dt',  function() { setTimeout(add_display_events,0);});
    $('#cnvTable').on('page.dt',  function() { setTimeout(add_display_events,0);});
    $('#cnvTable').on('search.dt',  function() { setTimeout(add_display_events,0);});
    
}


var highVariabilityRegions = [
  { chr: 'chr6', start: 27000000, end: 36600000 } // MHC / HLA regions
]


function showHVR(visible) {

    if(visible) {
        for(var i=0; i<highVariabilityRegions.length;++i) {
            var hvr = highVariabilityRegions[i];
            filters.push({expr:'!(chr=="'+hvr.chr+'" && start > ' + hvr.start + ' && end < ' + hvr.end + ') /* HVR */'}) 
        }
    }
    else {
        filters = $.grep(filters, function(f) {
                return f.id || f.expr.indexOf('HVR')<0;
        });
    }

    filterTable();
}

var highlightTr = null

function getCNVIndexForTableRow(tr) {
    return parseInt($(tr).find('a')[0].href.match(/cnv_([0-9]*)_detail/)[1],10)-1;
}

function showCNVDetailsForRow(tr) {
    var cnvIndex = getCNVIndexForTableRow(tr);
    console.log("Displaying CNV " + cnvIndex);
    show_cnv_details(cnvIndex);

    if(highlightTr)
       $(highlightTr).removeClass('highlighted');

    highlightTr = tr;
    $(highlightTr).addClass('highlighted');
}


function add_display_events() {

    console.log("Add display events");
    $("#cnvTable tbody tr").unbind('click').click(function(){
        console.log("Show cnv details");
        showCNVDetailsForRow(this);
     });
}

var TYPE_DESCRIPTIONS = {
    DEL : 'Deletion',
    DUP : 'Insertion',
    SNP : 'SNV'
}

function igv_load_cnv(cnv, pos) {
 var bam_url = 'http://localhost:60151/load?file='+encodeURIComponent(bam_file_path+bams[cnv.sample])+'&locus='+cnv.chr+'%3A'+cnv.start+'-'+cnv.end

 var pos_start
 var pos_end 
 if(pos) {
     pos_start = pos;
     pos_end = pos;
 }
 else {
     pos_start = cnv.start;
     pos_end = cnv.end;
 }
 
 $.ajax( bam_url, {
    complete: function() {
        igv_load(cnv.chr, pos_start, pos_end);
    }
 }); 
}

function igv_load(chr, start, end) {
    var locus_url = 'http://localhost:60151/goto?locus='+chr+'%3A'+start+'-'+end;
    $.ajax( bam_url, {
        complete: function() {
            $.ajax(locus_url);
        }
   });
}

function initAnnotations(cnvIndex) {
    if(!userAnnotations[cnvIndex])
        userAnnotations[cnvIndex] = { tags: { } };
}

var southLayout = null;

/**
 * Add 'chr' to a chromosome reference if necessary
 */
function ucscChr(chr) {
   return chr.startsWith('chr') ? chr : 'chr' + chr;
}

/** 
 * Show details of a CNV in the bottom pane
 */
function show_cnv_details(cnvIndex) {
    var cnv = cnvs[cnvIndex];
    var variants = cnvVariants[cnvIndex];
    var south = $('#ui-layout-south')

    if(cnvLayout.state.south.size<100)
        cnvLayout.sizePane("south",480);

    if(!southLayout) {
        south.html('<div class=ui-layout-center id=southwest></div><div class="ui-layout-east" id="southeast"></div>');
        southLayout = south.layout({ applyDefaultStyles: true, onresize: () => show_cnv_details(cnvIndex) });
        southLayout.sizePane("east",1024);
    }

    var southWest = $('#southwest');
    southWest.html('');
    
    var tagDiv = southWest.$div({id:"tags",'class':'cnvTags'});
    with(tagDiv) {
            
        var cnvTags = userAnnotations[cnvIndex] && userAnnotations[cnvIndex].tags;
        allTags.forEach(function(tag) {
            tagDiv.$div({'class':'cnvTag ' + ' tag'+(allTags.indexOf(tag)) + ' ' + (cnvTags && cnvTags[tag] ? 'assignedTag ' : 'unassignedTag')}).$span(tag).click(function() {
                confirm('Remove tag ' + tag + '?');
                delete userAnnotations[cnvIndex].tags[tag];
                addRowTags(highlightTr, userAnnotations[cnvIndex].tags);
                saveSettings();
                show_cnv_details(cnvIndex);
            });
        });
        
        $button('Tag').click(partial(addTagToRow,cnvIndex, show_cnv_details));
    }
    
    let sample = idMaskRegExp ? cnv.sample.match(idMaskRegExp)[1] : cnv.sample

    with(southWest.$table({id:"cnvDetails",'class':'cnvDetailTable'})) {
        with($tr()) { $th("Type"); $td(TYPE_DESCRIPTIONS[cnv.type]); }
        with($tr()) { $th("Sample"); $td(sample); }
        with($tr()) { $th("Genes"); $td(cnv.genes.join(",")); }
        with($tr()) { $th("Concordance"); $td(cnv.count); }
    }
    

    southWest.$p().$b("Region")
    with(southWest.$ul()) {
        $li().$a({href:'http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name='+cnv.chr+'%3A'+cnv.start+'-'+cnv.end+';search=Search', target:'ximmerdgv'}).$span('DGV')
        $li().$a({
                  href:'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+ucscChr(cnv.chr)+'%3A'+cnv.start+'-'+cnv.end+'&dgv=pack&knownGene=pack&omimGene=pack',
                  target: 'ximmerucsc'
                 }).$span('UCSC')
        with($li()) {
          with($a({href:'#'})) {
            $span('IGV')
            click(function() {
                igv_load_cnv(cnv);
                return false;
            });
          }
        }
    }

    southWest.$p().$b("Genes")
    var t = southWest.$table({id:'cnvTable'})
    $.each(cnv.genes, function() {
        var gene = this;
        with(t.$tr()) {
            $td(gene.toString())
            with($td()) {
                with($ul()) {
                    $li().$a({
                      href:'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+encodeURIComponent(gene)+'&search='+gene+'#diseases', 
                      target:'ximmergenecards'
                    }).$span('GeneCards');
                    $li().$a({href:'http://www.omim.org/search?index=entry&start=1&limit=10&search='+gene+'&sort=score+desc%2C+prefix_sort+desc', target:'ximmeromim'}).$span('OMIM');
                    if(geneList && geneList[gene])
                        $li().$span('Category ' + geneList[gene]);
                }
            }
        }
    })

    if(variants.length > 0) {
        with(southWest.$table({id:"cnvVariants",'class':'cnvDetailTable'})) {
            with($thead()) { with($tr()) { 
               $th({colspan:2, align:'center'}).$span("Variants");
               $th('Impact');
               $th('Severity');
               $th('Type');
               $th('Ref');
               $th('Alt');
            }}
            $.each(variants, function() {
                var v = this;
                with($tr()) {
                    var encPos = ucscChr(v.chr) +'%3A'+v.alleles[0].start;
                    with($td()) { 
                        $a({href:'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+encPos+'&dgv=pack&knownGene=pack&omimGene=pack', target:'ximmerucsc'}).$span(v.chr+':'+v.alleles[0].start)
                    }
                    with($td().$a({href:'#'})) { $span('igv'); click(function() { igv_load(v.chr, v.pos, v.alleles[0].start)}) }
                    $td(v.effect ? v.effect.toLowerCase() : 'Unknown')
                    $td(v.impact ? v.impact.toLowerCase() : 'Unknown');
                    $td(v.dosage == 1 ?  'Het' : 'Hom')
                    $td(v.depths[0])
                    $td(v.depths[1])
                }
            });
        }
    }
    
    var southEast = $('#southeast');
    southEast.html('');
    
    var plotWidth = Math.max(100,Math.round(southLayout.state.east.size-50));
    var plotHeight = Math.max(100,Math.round(southLayout.state.east.innerHeight-50));
    southEast.$div({id:'coveragePlot','class': 'coveragePlotDiv'});
    
    /*
    {
        $button({'id':'ampliconButton'}).$span('Amplicons');
        // $img({src:imgpath+'cnv_'+cnv.chr + '_'+cnv.start+'_'+cnv.end+'_' + cnv.sample + '.png', width:Math.max(100,Math.round(southLayout.state.east.size-50))})
       
        $svg({id:'plotSvg',
           width: plotWidth,
           height: plotHeight 
        });
        
        
    }
    */
    
    $('#coveragePlot').html(`<svg id=plotSvg width=${plotWidth} height=${plotHeight}></svg>`);
    
    var jsonPath = imgpath+'cnv_'+cnv.chr + '_'+cnv.start+'_'+cnv.end+'_' + cnv.sample + '.js';
    window.dx = new CNVDiagram($('#plotSvg')[0], jsonPath, plotWidth-50, plotHeight);
    dx.loadCnvData( () => {
        dx.renderPlot();
    });
    
    $('#ampliconButton').click(toggleAmplicons);
    window.cnvIndex = cnvIndex;
}

function addTagToRow(dataIndex, showFn) {
    var newTag = prompt('Enter the new tag: ','');
    if(!newTag)
        return;
    if(allTags.indexOf(newTag)<0)
        allTags.push(newTag);
                
    initAnnotations(dataIndex);
    userAnnotations[dataIndex].tags[newTag] = true;
    addRowTags(highlightTr, userAnnotations[dataIndex].tags);
    saveSettings();
    showFn(dataIndex);
}





