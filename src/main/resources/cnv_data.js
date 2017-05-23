/**
 * 
 */

            <% log.println "Writing cnv json" %>
            var cnvs = <%=groovy.json.JsonOutput.toJson(
                cnvs.collect { cnv ->
                  def variantInfo = [:];
                  if(hasVariants) {
                      variantInfo = [
                       het: cnv.het,
                       hom: cnv.hom, 
                       bal: ((cnv.bal == null || cnv.bal.isNaN())?null:cnv.bal), 
                       truncated: cnv.variants.any { it.maxEffect?.isTruncating() }
                      ]
                  };
                
                  return [ chr : cnv.chr, start: cnv.from, end: cnv.to ] + 
                      ['sample','type','count','stotal'].collectEntries { [it, cnv.getProperty(it)] } +
                      [ genes : cnv.genes ? cnv.genes.split(",") : [] ] +
                      cnv_callers.collect { caller -> 
                          if(cnv[caller])
                            [ caller, [ quality : cnv[caller].quality ] ] 
                          else
                            [ caller, false ]
                      }.collectEntries() + variantInfo
                }
            )%>;
            <% log.println "Writing variant json" %>
            var cnvVariants = <%="[" + cnvs.collect{ "["+ (it.variants*.toJson(it.sample)?:[]).join(",") + "]" }.join(",") + "]"%>;
            <% log.println "Writing callers json" %>
            var callers = <%=groovy.json.JsonOutput.toJson(cnv_callers)%>;
            <% log.println "Writing annotations json" %>
            <% def annotations = cnvs.collect { cnv -> 
                   cnvAnnotator ? cnvAnnotator.annotate(new Region(cnv.chr, cnv.from..cnv.to), anno_types[cnv.type]) : [spanning: 0f, spanningFreq: 0f] }; 
             %>
            var annotations = <%=groovy.json.JsonOutput.toJson(annotations)%>;
            
            var bam_file_path = <%=groovy.json.JsonOutput.toJson(bam_file_path)%>;

            var bams = <%=groovy.json.JsonOutput.toJson(sample_bams)%>;            
