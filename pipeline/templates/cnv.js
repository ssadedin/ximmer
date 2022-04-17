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

var initialized = false;

var cnvTable = null;
var cnvData = [];
var columns = [];

// Array of strings to 'eval' as filters
var filters = [];
var layout = null;

$(document).ready(function() {

    if(initialized)
        return;

    console.log("init");

    columns = []
    $('#cnvTable thead th').each(function() {columns.push(this.innerHTML);})
    columns = columns.map(function(x) { return {title: x}; })

    cnvTable = $('#cnvTable').DataTable({ 
        "iDisplayLength": 50,
        "columnDefs": [
            { "width": "2em", "targets": 0 },
            { "width": "25%", "targets": 0 },
            { "width": "6em", "targets": 0 },
            { "width": "4em", "targets": 0 }
        ]
    });

    cnvData = cnvTable.rows().data();

    layout = $('body').layout({ applyDefaultStyles: true });

    $(".cnvimg").click(function() {
        window.open(this.src, "cnvimg", "width=1600,height=1024");
    });

    $('#cnvTable').on('order.dt',  add_display_events);

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
        $button('Share').click(function() {
            alert('Copy the text below to share these filters with somebody else:\n\n'+JSON.stringify(JSON.parse(localStorage.cnvFilters)[location.href]));
        });

        if(typeof(localStorage.cnvFilters) != 'undefined') {
            $button('Import').click(function() {
                var newFilters = prompt("Paste filters from another session below:");
                if(newFilters) {
                    filters = JSON.parse(newFilters);
                    renderFilters();
                    filterTable();
                }
            });
        }
    }


    $('#addFilter').click(function() {
        filters.push({ expr: '', id: filters.length});
        renderFilters();
    });
    $('#runFilters').click(filterTable);

    var attrs = $.map(cnvs[0], function(v,k) { return k=='seqnames' ? 'chr' : k });
    $.each(annotations[9],function(k,v) { attrs.push(k); }); 

    $('#addFilter').one("click",function() {
        with($('#filterHelp')) {
            $span("Filter attributes: " + attrs.join(","));
        }
        layout.sizePane("north",150);
        $('#filterHelp').slideDown();
    });

    if(localStorage.cnvFilters) { 
        var oldFilters = JSON.parse(localStorage.cnvFilters)[location.href];
        if(oldFilters) {
            filters = oldFilters;
            renderFilters();
            filterTable();
        }
    }

    var cnvIndex = 1;
    $('h3').each(function() {
        $(this).$button('Filter Out').click(function() {
            var cnvId = this.parentNode.id.replace(/[a-z_]/g,'');
            console.log("filter out cnv " + cnvId);
            filters.push({expr: 'index!="'+cnvId+'"'}); 
            filterTable();
        });
    });

    initialized = true;
});

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

function updateFilters() {
    for(var i=0; i<filters.length; ++i) {
        if(filters[i].id != null)
            filters[i].expr = document.getElementById(filters[i].id).value;
    }
    var allFilters = localStorage.cnvFilters ? JSON.parse(localStorage.cnvFilters) : {};
    allFilters[location.href] = filters;
    localStorage.cnvFilters = JSON.stringify(allFilters);
}


var dataMap = {
    'TRUE' : true,
    'FALSE' : false
}

function filterTable() {

    updateFilters();

    // Update the filters from the fields
    // Filter out rows
    var newTable = cnvData.filter(function(value,index) { 
        var cnv = cnvs[index];
        cnv.chr = cnv.seqnames;
        cnv.index = index+1;
        window.dataval = value;
        var data = $.extend({}, cnv, annotations[index]);

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
    $('#cnvTable').dataTable({ data: newTable, 
                               columns: columns, 
                               iDisplayLength: 50,
                               destroy: true,
                               });
    add_display_events();
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

function add_display_events() {
    // $(".cnvRow").unbind('click').click(function() {
    $("#cnvTable tbody tr").unbind('click').click(function() {
        var cnvIndex = parseInt($(this).find('a')[0].href.match(/cnv_([0-9]*)_detail/)[1],10)-1;
        console.log("Displaying CNV " + cnvIndex);
        show_cnv_details(cnvIndex);
    })
}

var TYPE_DESCRIPTIONS = {
    DEL : 'Deletion',
    DUP : 'Insertion',
    SNP : 'SNV'
}

/** 
 * Show details of a CNV in the RHS pane 
 */
function show_cnv_details(cnvIndex) {
    var cnv = cnvs[cnvIndex];
    var variants = cnvVariants[cnvIndex];
    var east = $('#ui-layout-east')
    east.html('');

    with(east.$table({id:"cnvDetails",'class':'cnvTable'})) {
        with($tr()) {
            $th("Type");
            $td(TYPE_DESCRIPTIONS[cnv.type]) 
        }
        with($tr()) { $th("Genes"); $td(cnv.genes); }
        with($tr()) { $th("Concordance"); $td(cnv.count); }
    }

    with(east.$table({id:"cnvVariants",'class':'cnvTable'})) {
        with($thead()) { with($tr()) { $th({colspan:2, align:'center'}).$span("Variants") } }
        $.each(variants, function() {
            var v = this;
            with($tr()) {
                with($td()) { 
                    $a({href:'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+v.chr+'%3A'+v.alleles[0].start+'&dgv=pack&knownGene=pack&omimGene=pack'}).$span(v.chr+':'+v.alleles[0].start)
                    $span(' ')
                    $a({href:'http://localhost:60151/goto?locus='+v.chr+'%3A'+v.alleles[0].start}).$span('igv');
                }   
            }
        });
    }
    east.$p().$b("Region")
    with(east.$ul()) {
        $li().$a({href:'http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19/?name='+cnv.seqnames+'%3A'+cnv.start+'-'+cnv.end+';search=Search'}).$span('DGV')
        $li().$a({href:'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position='+cnv.seqnames+'%3A'+cnv.start+'-'+cnv.end+'&dgv=pack&knownGene=pack&omimGene=pack'}).$span('UCSC')
        $li().$a({href:'http://localhost:60151/goto?locus='+cnv.seqnames+'%3A'+cnv.start+'-'+cnv.end}).$span('IGV')
    }

    east.$p().$b("Genes")
    var t = east.$table({id:'cnvTable'})
    $.each(cnv.genes.split(","), function() {
        var gene = this;
        with(t.$tr()) {
            $td(gene.toString())
            with($td()) {
                with($ul()) {
                    $li().$a({href:'http://www.genecards.org/cgi-bin/carddisp.pl?gene='+encodeURIComponent(gene)+'&search='+gene}).$span('GeneCards');
                    $li().$a({href:'http://www.omim.org/search?index=entry&start=1&limit=10&search='+gene+'&sort=score+desc%2C+prefix_sort+desc'}).$span('OMIM');
                }
            }
        }
    })
}

