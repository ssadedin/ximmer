
printf = function(x,...) { print(sprintf(x,...)) }

sjia.analysis.names.replaced = list(
  "comparesim0.1" = "exome_ex1",
  "comparesim0.2" = "exome_ex1",
  "comparesim3.1" = "comparesim1",
  "comparesim3.2" = "comparesim3.2",
  "comparesim5.1" = "exome_ex5",
  "comparesim5.2" = "exome_ex5"
)

# correct, downsampled reads
sjia.analysis.names.ds = list(
  "comparesim0.1" = "exome_ex1.ds",
  "comparesim0.2" = "exome_ex1.ds",
  "comparesim3.1" = "comparesim1ds",
  "comparesim3.2" = "comparesim3.2ds",
  "comparesim5.1" = "exome_ex5.ds",
  "comparesim5.2" = "exome_ex5.ds"
)

mg.analysis.names.untrimmed= list(
  "sim1.1" = "mgex",
  "sim1.2" = "sim1.2",
  "sim3.1" = "sim3.1",
  "sim3.2" = "sim3.2",
  "sim5.1" = "sim5.1",
  "sim5.2" = "sim5.2"
)

# Small analysis based on trimmed data
# Note this data does not include Excavator results
# ie: sim.callers = c("ed","xhmm","mops")
mg.analysis.names.trimmed = list(
  "sim1.1" = "sim1.1",
  "sim3.1" = "sim3.1",
  "sim5.1" = "sim5.1"
)

mg.analysis.names.final = list(
  "sim1.1" = "sim1.1",
  "sim1.2" = "sim1.2",
  "sim3.1" = "sim3.1",
  "sim3.2" = "sim3.2",
  "sim5.1" = "sim5.1",
  "sim5.2" = "sim5.2"
)



halo_sim_loaders = list(
  ex=function(batch,sims) { do.call(c,lapply(sims,function(sim) load_excavator_results(sprintf("%d/sim%d.%d/source_files/sim%d.%d.excavator.cnvs.tsv",sim$size,sim$size,sim$rep,sim$size,sim$rep),sprintf("%s_%s",sim$size,sim$rep))))},
  ed=function(batch,sims) {do.call(c,lapply(sims,function(sim) load_exomedepth_results(sprintf("%d/sim%d.%d/source_files/sim%d.%d.exome_depth.cnvs.tsv", sim$size, sim$size,sim$rep,sim$size,sim$rep),sprintf("%s_%s",sim$size,sim$rep))))},
  xhmm=function(batch,sims) {do.call(c,lapply(sims,function(sim) load_xhmm_results(Sys.glob(sprintf("%d/sim%d.%d/source_files/*.xhmm_discover.xcnv",sim$size,sim$size,sim$rep)[[1]]),sprintf("%s_%s",sim$size,sim$rep))))},
  cnmops=function(batch,sims) {do.call(c,lapply(sims,function(sim) load_cn_mops_results(Sys.glob(sprintf("%d/sim%d.%d/source_files/*.cnmops.cnvs.tsv",sim$size,sim$size,sim$rep)[[1]]),sprintf("%s_%s",sim$size,sim$rep))))},
  cfr=function(batch,sims) {do.call(c,lapply(sims,function(sim) load_conifer_results(Sys.glob(sprintf("%d/sim%d.%d/source_files/*.conifer.cnvs.tsv",sim$size,sim$size,sim$rep)[[1]]),sprintf("%s_%s",sim$size,sim$rep))))},
  anghmm=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) {
      pattern = sprintf("%d/sim%d.%d/%s/sim%d.%d.*counts.angelhmm*.cnvs.bed", sim$size,sim$size,sim$rep,batch,sim$size,sim$rep)
      print(sprintf("Pattern = %s", pattern))
      load_angel_results(
        Sys.glob(pattern)[[1]],sprintf("%s_%s",sim$size,sim$rep))
    }))}
)

mg.analysis.names.ft = mg.analysis.names.trimmed

sim.caller.loaders = list(
  ex=function(batch,an, sims,num) { do.call(c,lapply(sims,function(sim) load_excavator_results(sprintf("%s/%s.%d/excavator/%s.%d.excavator.cnvs.tsv",sim,an[[sim]],num,an[[sim]],num))))},
  ec=function(batch,an, sims,num) {do.call(c,lapply(sims,function(sim) load_exome_copy_results(sprintf("%s/%s.%d/exome_copy/%s.%d.exome_copy.cnvs.tsv",sim, an[[sim]],num,an[[sim]],num))))},
  ed=function(batch,an, sims,num) {do.call(c,lapply(sims,function(sim) load_exomedepth_results(sprintf("%s/%s.%d/exome_depth/%s.exome_depth.cnvs.tsv",sim, an[[sim]],num,an[[sim]]))))},
  xhmm=function(batch,an, sims,num) {do.call(c,lapply(sims,function(sim) load_xhmm_results(sprintf("%s/%s.%d/xhmm/%s.%d.params.xhmm_discover.xcnv",sim, an[[sim]],num,an[[sim]],num))))},
  mops=function(batch,an, sims,num) {do.call(c,lapply(sims,function(sim) { load_cn_mops_results(Sys.glob(sprintf("%s/%s.%d/cn_mops/*.cn_mops_call_cnvs.tsv",sim, an[[sim]],num,an[[sim]]))[[1]])}))}
)

compute.sim.info = function(sim.all, batch.name) {
   lapply(sim.all, function(sim) {
    sim$report=sprintf("%d/sim%d.%d/%s/sim%d.%d.cnv.chrX.summary.tsv", sim$size, sim$size,sim$rep, batch.name, sim$size,sim$rep)
    sim$name=sprintf("%d_%d",sim$size,sim$rep)
    return(sim)
  }) 
}

load.halo.combined.results = function(sims, callers=sim.callers) {
  result.table = do.call(rbind, unname(lapply(sims, function(sim) {
    # Load the results
    results.name = sim$report
    printf("Loading from file %s", results.name)
    r = read.table(results.name,header=T)
    r$sample= paste(r$sample,sprintf("%d_%d",sim$size,sim$rep),sep="_")
    r$sim = sim$name
    return(r)
  })))
  
  results = GRanges(seqnames=result.table$seqnames, ranges=IRanges(result.table$start, result.table$end),
                    truth=result.table$truth,
                    sample=result.table$sample,
                    sim=result.table$sim)  
  
  for(caller in callers) {
    mcols(results)[,caller] = result.table[,caller]
  }
  return(results)
}

load.combined.results = function(sims, batch.name, pattern="%s/%s/report/%s.cnv.summary.tsv") {
  sim.names = names(sims)
  result.table = do.call(rbind, unname(lapply(sim.names, function(sim.name) {
    # Load the results
    results.name = sprintf(pattern,sim.name,batch.name, batch.name)
    printf("Loading from file %s", results.name)
    r = read.table(results.name,header=T,sep='\t')#[,1:15]
    r$sim =sim.name
    return(r)
  })))
  
  results = GRanges(seqnames=result.table$chr, ranges=IRanges(result.table$start, result.table$end),
                    truth=result.table$truth,
                    sample=result.table$sample,
                    sim=result.table$sim)  
  for(caller in sim.callers) {
    # workaround for egregious problem with cnmops being called mops some places
    result.caller = ifelse(caller == "mops", "cnmops", caller)
    mcols(results)[,caller] = result.table[,result.caller]
  }
  return(results)
}

read_cnv_params = function(file.name) {
  params.raw = strsplit(read.table(file.name,stringsAsFactors = F)$V2,"=")
  params = lapply(params.raw, function(p) p[[2]])
  names(params) = lapply(params.raw, function(p) p[[1]])  
  return(params)
}

get_param_labels = function(file.name) {
  printf("Reading %s", file.name)
  params = read_cnv_params(file.name)
  return(list(
    ex=paste0("theta=",params$excavator_theta),
    ec=sprintf("Q%s",params$exome_copy_quality_threshold),
    ed=sprintf("trans. prob=%s",params$transition_probability),
    xhmm=sprintf("cnv rate=%s",params$exome_wide_cnv_rate),
    mops=sprintf("prior impact=%s",params$prior_impact)
  ))
}

load.exome.result.set = function(params, truth, analysis.names, exclude.bed=NULL) {
  # Load a set of results based on the name of the run (eg: "optimal", or "default")
  # Looks up the correct run to use based on the values found in the params file 
  # named by the 'params' parameter. The results for this run are then used. This is
  # done for all callers in 'sim.callers' and a list is returned.
  results_for_params = list()
  sim.names = names(analysis.names)
  for(s in sim.callers) {
    # find the optimal result set for this caller
    for(i in 1:6) {
      if(get_param_labels(sprintf("params.%d.txt",i))[[s]]==params[[s]]) {
        printf("Found requested parameters in params.%d.txt for caller %s", i,s)
        results_for_params[[s]] = load_ranked_run_results(truth, analysis.names, exclude.bed, sim.names, i)[[s]]
      }
    }
  }
  return(results_for_params)
}

load_ranked_run_results_for_simulation = function(analysis.names) {
  # Loads the given set of simulation runs that are specified by 
  # analysis.names.
  #
  # simulation runs in the result set. The results are specified by
  # the analysis.names list. The names of the list specify top level
  # directories (specifyign simulations) of the simulation hierarchy 
  # and the values are the  actual run batch within that directory to load 
  # for each simulation
  sim.names = names(analysis.names)
  truth = load.truth(sim.names)
  results.all = lapply(1:6, function(j) {
    load_ranked_run_results(truth, analysis.names, exclude.bed=dgv.exclude, sim.names,j)
  })
}

find.tp.overlaps.opt = function(cnvs, truth) {
  cnvs$tpid = sapply(cnvs, function(cnv) {
    sample.tps = truth.by.sample[[cnv$sample]] 
    tps = findOverlaps(cnv,sample.tps)
    if(length(tps)>0) {
      result = sample.tps$index[subjectHits(tps)][[1]]
      return(result)
    }
    return(NA)
  })
  return(cnvs)
}

find.tp.overlaps.opt3 = function(cnv.calls, truth) {
  if(length(cnv.calls)==0)
    return(cnv.calls)
  
  cnv.calls$index = 1:length(cnv.calls)
  cnv.calls.df = as.data.frame(cnv.calls)
  truth.df = as.data.frame(truth)
  overlaps = findOverlaps(cnv.calls,truth)
  cnv.calls$tpid = sapply(cnv.calls, function(cnv) {
    tp.overlaps = overlaps[queryHits(overlaps)==cnv$index & truth.df[subjectHits(overlaps),]$id==cnv$sample]
    tps = truth.df[subjectHits(tp.overlaps),]
    if(nrow(tps)>0) {
      result = tps[1,]$index
      return(result)
    }
    return(NA)
  })
  return(cnv.calls)
}

load_ranked_run_results = function(truth, analysis.names, exclude.bed=NULL, batch=batch.name, callers=sim.callers, filter.samples=NULL, filterChrX=T, ...) {
  
  if(!is.null(filter.samples)) {
    truth = truth[filter.samples(truth$id)]
  }
  
  # Start by setting the index of each true positive (we need it later)
  truth$index = 1:length(truth)
  
  sim.results = lapply(callers, function(caller) {
    
    if(is.null(sim.caller.loaders[[caller]]))
      stop(sprintf("A CNV caller %s was specified in the callers argument but is not found in the list of loaders: %s",
                   caller, paste0(names(sim.caller.loaders),collapse=",")))
    
    cnvs = sim.caller.loaders[[caller]](batch,analysis.names,...)
    
    if(!is.null(filter.samples)) {
      cnvs = cnvs[filter.samples(cnvs$sample)]
    }
    
    printf("Loaded %d raw CNVs for %s", length(cnvs), caller)
    
    if(filterChrX) {
      cnvs = cnvs[cnvs$type == "DEL" & seqnames(cnvs)=='chrX']
      seqlevels(cnvs) = "chrX"
    }
    
    printf("%d deletion CNVs for %s", length(cnvs), caller)
    
    # Remove any false positives that are over the excluded region
    # (if it is provided)
    if(!is.null(exclude.bed)) {
      cnvs = cnvs[ !(cnvs %over% exclude.bed) | (cnvs %over% truth)  ]
    }
    
    cnvs$sample = as.character(cnvs$sample)
    cnvs$true = unlist(lapply(cnvs, function(cnv) {
      cnv %over% truth[truth$id == cnv$sample]
    }))
    
    # What if a caller produces multiple calls over the same true positive?
    # We should not count them as multiple TP calls!
    # We assign a "true positive id", tpid to each true positive. Then 
    # we can use that later to ensure we don't count a true positive more
    # than once
    printf("Finding true positives overlapping %d called cnvs for %s",length(cnvs), caller)
    x <<- cnvs
    cnvs = find.tp.overlaps.opt3(cnvs,truth)
    
    # order by quality so that we count true positives correctly to produce ROC-style curve
    printf("Counting true positives ordered by rank for %s",caller)
    cnvs = cnvs[order(cnvs$qual,decreasing = T)]
    
    if(length(cnvs)>0) {
      true.count = 0
      true.counts = c()
      tp.identified = c()
      call.count = 0
      calls = c()
      
      for(i in 1:length(cnvs)) {
        print(sprintf("CNV %d by caller %s", i, caller))
        cnv = cnvs[i]
        if(cnv$true && (cnv$tpid %in% tp.identified)) {
          printf("CNV %d identified multiple times by different calls",i)        
        }
        else
        if(cnv$true) {
            true.count = true.count + 1
            call.count = call.count + 1
            tp.identified = c(tp.identified, cnv$tpid)
        }
        else  
          call.count = call.count + 1
        
        calls = c(calls, call.count)
        true.counts = c(true.counts, true.count)
      }
      cnvs$count = true.counts
      cnvs$index = 1:length(cnvs)
      cnvs$calls = calls
    }
    else {
      cnvs$count = c()
      cnvs$index = c()
      cnvs$calls = c()
    }
    return(cnvs)
  })
  names(sim.results)=callers
  sim.results$truth = truth
  return(sim.results)
}

plot_caller_roc = function(truth, caller.result, numbers=F,x.offset=0,...) {
  # x = 1-(caller.result$count/caller.result$index)
  fp = (caller.result$calls-caller.result$count)
  # x = fp / (length(caller.result) - max(caller.result$count))
  x = fp # / caller.result$index
  
  y = caller.result$count # /length(truth)
  # x.plot = c(0,jitter(x))
  x.plot = c(0,x) + x.offset
  # y.plot = c(0,jitter(y))
  y.plot = c(0,y)
  
  lines(x.plot,y.plot,...)
  # points(x.plot,y.plot,pch=19,cex=0.2,...)
  if(numbers)
    text(x,y,paste0(caller.result$count,"/",caller.result$index))
}

roc_plot= function(legend.caller, position, title, results, truth, param.x.offset=0.025,figletter=F, ...) {
  max.x = max(c(unlist(sapply(results, function(result) result[[legend.caller]]$index - result[[legend.caller]]$count))),5)
  # max.y = max(c(unlist(sapply(results, function(result) result[[legend.caller]]$count))),5)
  max.y = length(truth)
  plot(c(0), xlim=c(0,max.x), ylim=c(0,max.y),type="n", xlab="False Positives", ylab=sprintf("True Positives found out of %d Total",max.y))
  params = lapply(1:6,function(j) get_param_labels(sprintf("params.%d.txt",j))[[legend.caller]])
  
  plot.order = roc_plot_order(results,legend.caller)  
  legend(position,fill=plot.order, legend=params[plot.order])
  title(main=title)
  x_per_inch = max.x / par()$pin[1]
  y_per_inch = max.y / par()$pin[2]
  
  plot.offset= max(-1, -length(plot.order)*param.x.offset*x_per_inch/2)
  for(j in plot.order) {  
    plot_caller_roc(truth, results[[j]][[legend.caller]], col=j, F, plot.offset,...)  
    plot.offset = plot.offset+ param.x.offset*x_per_inch
  }
  mtext(sprintf("True Positives found from %d True Deletions", length(truth)),line=0.8)  
  if(figletter!=F) {
    text(x_per_inch*0.3,max.y - y_per_inch*0.5 ,figletter,cex=2, font=2)
  }
}

roc_plot_order = function(results, caller) {
  # It looks better on the plot if the callers / configs with less true positives 
  # appear on top of the ones with more, so this function returns an ordering
  # to achieve that
  order(unlist(sapply(results, function(r) max(r[[caller]]$count))),decreasing = F) 
}

load.truth = function(sims) {
  # Load the "true_cnvs.bed" file from each simulation which is expected to 
  # be in the directories named by the "sims" list given as argument.
  truth = do.call(c,lapply(sims, function(s) {
    sim.truth = read.bed.ranges(sprintf("%s/true_cnvs.bed",s))
    sim.truth$source = s
    return(sim.truth)
  }))
  truth$id = as.character(truth$id)  
  return(truth)
}

plot.exome.result.set = function(optimal.results,
                                 truth,
                                 newPlot=T, 
                                 plot.offset=0,
                                 plot.offset.increment=0,
                                 cols=nice_colors,
                                 plot.lwd.increment=0,
                                 plot.lwd=3,
                                 legend.cex=1.0,
                                 ...) {
  
  max.x = max(c(unlist(sapply(optimal.results, function(result) result$index[[length(result)]] - result$count[[length(result)]]))),5)
  max.y = length(truth)
  x_per_inch = max.x / par()$pin[1]
  
  if(plot.offset.increment == 0)
    plot.offset.increment = 0.06*x_per_inch
  
  if(newPlot) {
    plot(c(0), xlim=c(0,max.x), ylim=c(0,max.y),type="n", xlab="False Positives", ylab=sprintf("True Positives Detected (out of %d total)",max.y),...)
  }
  for(s in sim.callers) {
    plot_caller_roc(truth, optimal.results[[s]], col=cols[[match(s,sim.callers)]], F, x.offset=plot.offset, lwd=plot.lwd, ...)
    plot.offset = plot.offset+ plot.offset.increment
    plot.lwd = plot.lwd + plot.lwd.increment
    printf("Plot offset = %f, Width=%f", plot.offset, plot.lwd)
  }
  legend("topright", fill=cols[1:length(sim.callers)], legend=lapply(sim.callers,function(caller) sim.caller.labels[[caller]]), cex=legend.cex)
}

compute_deletion_metrics = function(deletion.ranges, exome_bed_file) {  # "../design/EXOME.bed"
 
  if(class(deletion.ranges)[[1]]!="GRanges")
    stop("Deletion ranges should be a GRanges object")
 
  # Load the exomes
  exome = reduce(read.bed.ranges(exome_bed_file))
  
  # For each simulated CNV, calculate number of target regions included
  print("Computing targets")
  deletion.ranges$targets = unlist(unname(lapply(deletion.ranges, function(r) {
    length(unique(subjectHits(findOverlaps(r,exome))))
  })))
  
  print("Computing seqbp")
  deletion.ranges$seqbp = unlist(unname(lapply(deletion.ranges, function(r) {
    targets = exome[unique(subjectHits(findOverlaps(r,exome)))]
    sum(width(targets))
  })))
  
  print("Computing bp")
  deletion.ranges$bp = unlist(unname(lapply(deletion.ranges, function(r) {
    targets = exome[unique(subjectHits(findOverlaps(r,exome)))]
    max(end(targets)-min(start(targets)))
  })))  
  return(deletion.ranges)
}


calc.binned.result.perf = function(analysis.names, results.metrics, exome.bed, bin.by, bin.levels, exclude=NULL, calc) {
  #
  # Load all the simulation results into one giant table
  #
  # results -    The results object which should be a list of GRanges objects,
  #              named by caller. These are loaded by "load_ranked_run_results"
  #              functions
  # exome.bed  - string specifying location of bed file of exome target regions
  # bin.by     - string specifying a valid column on the results object that will
  #              be loaded to bin by. Columns include anything natively on the results,
  #              (loded by load.combined.results), but also anything added by
  #              compute_deletion_metrics, including "targets" and "seqbp"
  # bin.levels - the set of threshold to bin using. Must cover range of bin.by.
  # 
  
  lapply(results.metrics, function(results.raw) {
    
    # Compute the size metrics of all the deletions in the results
    printf("Computing deletion metrics for %d results", length(results))
    printf("Calculating performance for %d bins", length(bin.levels))
    
    # if exclusion was provided, remove all results that overlap the exclusion
    results = results.raw[ !(results.raw %over% exclude) ]
    
    bin.by.values = mcols(results)[,bin.by]
    #printf("Unique binning values are: %s", capture.output(print(unique(bin.by.values))))
    
    if(max(bin.by.values) > max(bin.levels))
      stop(sprintf("The range of bin.levels does not fully encompass the range of the values in bin.by (max=%d): '%s'", max(bin.by.values), bin.by))
    
    results.bins = cut(bin.by.values, bin.levels)
    
    # r = results[results.bins == levels(results.bins)[1]]
    
    results.binned = by(as.data.frame(results), results.bins, function(r) {
      calc(r)
    })
    #bin.counts = by(as.data.frame(results), results.bins, function(r) {
    #  length(r)
    #})
    #print(bin.counts)
    #printf("In bins %s ", str(bin.levels), " counts are %s", str(bin.counts))
  })
}

sens_by_del_plot_frame = function(bin.levels=F, 
                                  log.scale=F, 
                                  plot.title="Sensitivity vs Simulated Deletion Size", 
                                  plot.xlab="Deletion Size (sequenced bp)", 
                                  plot.ylab="Sensitivity (fraction of true CNVs discovered)",
                                  x.max=0,
                                  cex.lab=1.0,
                                  draw.axis=T
) {
  binned.perf.plot.frame(bin.levels, log.scale, plot.title, plot.xlab, plot.ylab, x.max, cex.lab, draw.axis)
}

binned.perf.plot.frame = function(bin.levels=F, log.scale=F, plot.title, plot.xlab, plot.ylab, x.max=0, cex.lab=1.0, draw.axis=T) {
  
  if(log.scale) 
    scale.fn = log10
  else
    scale.fn = function(x) { x }
  
  if(x.max == 0) {
    x.max = scale.fn(max(as.numeric(bin.levels)))
    printf("Scaling x.max automatically using highest bin")
  }
  
  printf("X.max = %f", x.max)
  plot(0, xlim=c(0,x.max), 
       ylim=c(0,1), main=plot.title, 
       xlab=plot.xlab, ylab=plot.ylab,
       type="n",
       xaxt="n",
       cex.lab=cex.lab
  )
  
  abline(h=(0:5)/5, col="#bbbbbb", lty=2); abline(v=(0:8) * round(x.max/8), col="#bbbbbb", lty=2) 
  
  if(draw.axis)
    axis(1,at=scale.fn(c(bin.levels[-1],x.max)), labels=c(bin.levels[-1],x.max))
}

plot.binned.performance = function(r,bin.levels,pch=19,points.cex=1,palette=nice_colors,log.scale=F, center.bins=T, 
                                   legend.pos="bottomright", 
                                   callers=sim.callers,
                                   caller.labels=sim.caller.labels,
                                   ...) {
  # 
  # Plots a line graph showing sensitivity vs bin
  #
  #   bin.levels - series of thresholds of bin boundaries, not including 0
  #   sim.callers - MUST BE SET AS GLOBAL VARIABLE
  #   r - list of lists of sensitivity, with names being the callers
  #
  if(log.scale) 
    scale.fn = log10
  else
    scale.fn = function(x) { x }
  
  for(i in 1:length(callers)) {
    print(callers[[i]])
    sens = c(0, unname(r[[callers[[i]]]]))
    
    if(center.bins) {
      levels.x = c(0, scale.fn(bin.levels[1:length(bin.levels)-1] + diff(bin.levels)/2))
    }
    else {
      levels.x = bin.levels
    }
    
    printf("bins = %s", capture.output(print(bin.levels)))
    printf("Perf = %s", capture.output(print(sens)))
    
    sens = sens[1:length(levels.x)]
    
    lines(levels.x, sens, 
          col=palette[[i]],
          lwd=0.8, 
          ...)
    points(levels.x, sens, 
           pch=pch,
           cex=points.cex,
           col=palette[[i]]
    )
  }
  
  # Legend
  legend(legend.pos,
         legend= paste(sapply(callers,function(s) caller.labels[[s]])),  
         pch=c(rep(19,length(callers)), rep(18, length(callers))), 
         col=rep(palette[1:4],18), 
         bg="white")
}


plot.binned.sensitity = plot.binned.performance



count.by.size.bin = function(results, truth.metrics, bin.levels) { 
  # Segregate the given true deletions in truth.metrics by sequenced base pairs
  # (deletion size) and return a GRanges object for each bin indicating for 
  # every caller whether the caller found the true positives in that size range.
  # 
  # The result object must have a 'sim' column that specifies which simulation or
  # other result set the result belongs to. This is compared to the truth.metrics$source
  # column. These must match for the two to be considered a matching TP.
  #
  
  truth.bins = cut(truth.metrics$seqbp, bin.levels)
  
  bin.cnv.results(truth.metrics, truth.bins, results)
#   
#   binning.result = by(as.data.frame(truth.metrics), truth.bins, function(truth.bin) {
#     #print(truth.bin)
#     truth.bin.results = mapply(function(cnv.chr, cnv.start, cnv.end, cnv.id, sim) {
#         sample.results = results[results$sample==cnv.id & results$sim==sim]
#         print(sprintf("Number of sample results = %d", length(sample.results)))
#         cnv.results = subjectHits(findOverlaps(GRanges(cnv.chr, IRanges(cnv.start,cnv.end)), sample.results))
#         if(length(cnv.results)>0) {
#           return(sample.results[cnv.results[[1]]])
#         } else {
#           return(GRanges(cnv.chr,IRanges(cnv.start,cnv.end), ed=FALSE, ex=FALSE, xhmm=FALSE, mops=FALSE, truth=TRUE, sample=cnv.id))        
#         }
#     }, truth.bin$seqnames, truth.bin$start, truth.bin$end, truth.bin$id, truth.bin$source)
#   })
#   
#   x<<-binning.result
#   
#   # Combine separate GRanges objects to 1 for each bin
#   return(lapply(binning.result, function(r.bin) {
#     do.call(c, unname(r.bin))
#   }))
}

bin.cnv.results = function(truth, truth.bins, results) { 
  #
  # Bins the given CNVs in 'truth' by the given bins and returns the given results partitioned into those
  # bins
  #
  binned = by(as.data.frame(truth), truth.bins, function(truth.bin) {
    
    results.cols = names(mcols(results))
    
    # print(truth.bin)
    truth.bin.results = mapply(function(cnv.chr, cnv.start, cnv.end, cnv.id, sim) {
      sample.results = results[results$sample==cnv.id & results$sim == sim]
      cnv.results = subjectHits(findOverlaps(GRanges(cnv.chr, IRanges(cnv.start,cnv.end)), sample.results))
      if(length(cnv.results)>0) {
        return(sample.results[cnv.results[[1]]])
      } else {
        result = GRanges(cnv.chr,IRanges(cnv.start,cnv.end))
        for(i in results.cols) {
          mcols(result)[,i] = F
        }
        # ed=FALSE,ex=FALSE,xhmm=FALSE,mops=FALSE,cfr=FALSE,truth=TRUE,
        result$sample=cnv.id
        result$sim=sim
        return(result);
      }
    }, truth.bin$seqnames, truth.bin$start, truth.bin$end, truth.bin$id,truth.bin$source)
  })
  
  # The above produces bins with a list of GRanges in each bin.
  # It's nicer to coalesce each list into a single GRanges with multiple ranges
  # Combine separate GRanges objects to 1 for each bin  
  binned.combined = lapply(binned, function(r.bin) { 
    if(length(r.bin)>0)
      do.call(c, unname(r.bin)) 
    else 
      return(r.bin)      
  })
}



 
plot.qscore.calibration = function(caller, qscore.results.all, newPlot=T, plot.title=NA, cex.lab=1.0, ...) {
  
  qscore.results = qscore.results.all[[caller]]
  
  phred.scaled.callers = c("xhmm","ed")
  
  if(caller %in% phred.scaled.callers) {
    qscore.results.bins = c(-10,0,5,10,20,30,40,50,60)
    plot.xlab = "Phred Scaled Quality Score"
  }
  else {
    qscore.results.bins = pretty(qscore.results$qual, 8)
    plot.xlab = "Raw Quality Score"
  }
  
  if(is.na(plot.title)) {
   plot.title = sprintf("Quality Score Calibration for %s", sim.caller.labels[[caller]])
  }
  
  qscore.results.bins.tprate = by(mcols(qscore.results), cut(qscore.results$qual, qscore.results.bins), function(cnvs) {
    printf("sum(cnvs$true)=%d",sum(cnvs$true))
    printf("length(cnvs)=%d",nrow(cnvs))
    return(sum(cnvs$true) / nrow(cnvs))
  })
  
  qscore.results.bins.tprate[is.na(qscore.results.bins.tprate)] = 0
  
  bin.centres = qscore.results.bins[-1]-(diff(qscore.results.bins)/2)
  
  if(newPlot) {
    plot(c(-10,bin.centres),
         c(0, qscore.results.bins.tprate),
         xlim=c(min(qscore.results.bins),max(qscore.results.bins)),
         ylim=c(0,1),
         pch=19,
         main=plot.title,
         xlab=plot.xlab,
         ylab="Proportion of true positives (Precision)",
         cex.lab=cex.lab,
         type="n"
    )
    grid()    
  }
  points(c(-10,bin.centres), c(0, qscore.results.bins.tprate), pch=19, ...)  
  lines(c(-10,bin.centres), c(0, qscore.results.bins.tprate), ...)
}



plot.halo.result.set = function(
  ranked.results,
  callers=sim.callers,
  newPlot=T, 
  plot.offset=0,
  plot.offset.increment=-1,
  cols=nice_colors,max.x=0, 
  plot.title="Deletion Detection Performance\nTrue Positives vs False Positives Ranked by Quality",
  title.cex=1.0,
  lty=1,
  legend=T,
  legend.cex=1.0,
  truth.set = NULL,
  ...) {
  
  if(is.null(truth.set)) {
    if(!is.null(ranked.results$truth)) { # if truth set attached to results, use that
      truth.set = ranked.results$truth
    }
    else {
      truth.set = truth # use global truth variable
    }
  }
  
  if(max.x==0)
    max.x = max(c(unlist(sapply(ranked.results, function(result) result$calls[[length(result)]] - result$count[[length(result)]]))),5)
  
  max.y <- length(truth.set)
  tmp.max.y = length(truth.set)
  
  print(sprintf("The y axis max = %d", length(truth.set)))
  print(sprintf("The y axis max = %d", max.y))
  print(sprintf("The tmp y axis max = %d", tmp.max.y))
  
  x_per_inch = max.x / par()$pin[1]
  
  if(plot.offset.increment < 0)
    plot.offset.increment = 0.05*x_per_inch
  
  if(newPlot) {
    plot(c(0), xlim=c(0,max.x), ylim=c(0,tmp.max.y),type="n", 
         xlab="False Positives", 
         ylab=sprintf("True Positives Found (out of %d Total)", tmp.max.y), 
         ...)
  }
  
  if(length(callers)==0) {
    stop("No CNV callers were listed in the callers argument")
  }
  
  for(i in 1:length(callers)) {
    s = callers[[i]]
    caller.lty = lty[[((i-1) %% length(lty))+1]]
    print(sprintf("caller=%s, lty = %d", s, caller.lty))
    
    if(is.null(ranked.results[[s]]))
      stop(sprintf("A CNV caller %s was specified for plotting, but no results were found in the ranked.results parameter for that name", s))
    
    col.index = match(s,callers)
    if(col.index>length(cols))
      stop(sprintf("There are more CNV callers defined (at least %d) than colors in the given color list (%s)", 
                   col.index, paste(cols,collapse = ",")))
    
    plot_caller_roc(truth.set, ranked.results[[s]], col=cols[[col.index]], F, x.offset=plot.offset, lwd=3, lty=caller.lty, ...)
    plot.offset = plot.offset+ plot.offset.increment
    printf("Plot offset = %f", plot.offset)
  }
  
  legend.labels = lapply(callers, function(caller) {
    if(caller %in% names(sim.caller.labels))
      sim.caller.labels[[caller]]
    else
      caller
  }) 
  
  if(legend) {
    legend("bottomright", col=cols[1:length(callers)], lty=lty, lwd=4, cex=legend.cex, legend=legend.labels)
  }
  title(main=plot.title, cex.main=title.cex)
}


halo.simulation.load.truth = function(sim.all) {
   truth = do.call(c,lapply(sim.all, function(s) {
    raw = read.bed.ranges(sprintf("%d/sim%d.%d/source_files/true_cnvs.bed",s$size,s$size,s$rep))
    raw$id= paste(raw$id,sprintf("%d_%d",s$size,s$rep),sep="_")
    # raw$source= paste(raw$id,sprintf("%d_%d",s$size,s$rep),sep="_")
    raw$source= sprintf("%d_%d",s$size,s$rep)
    return(raw)
  }))
  truth$id = as.character(truth$id)
  return(truth)
}

load.binned.cnv.truth.set = function(analysis.names, 
                                     batch.name, 
                                     truth, 
                                     bin.by, 
                                     bin.levels, 
                                     calc=cnv.caller.recall, 
                                     callers=sim.callers,
                                     filter.samples=NULL
                                     ) {
  
  # Load all the simulation results into one giant table
  # and bin them by cnv size
  #
  # Note: if analysis.names is not a character type then it will be treatedas the actual results 
  # as a GRanges object, as would be loaded by load.combined.results
  
  if(typeof(analysis.names)=="character") {
    print("Loading results from file ...")
    results = load.combined.results(analysis.names, batch.name)
  }
  else {
    print("Results are pre-loaded, using object directly")
    results = analysis.names
  }
  
  print("Cutting bins ...")
  truth.bins = cut(bin.by, bin.levels)
  
  print("Binning true CNVs ...")
  results.binned = bin.cnv.results(truth, truth.bins, results)
  
  result.tmp <<- results.binned
  
  # Compute sensitivity for each bin
  result.binned.sensitivity = lapply(callers, function(sim.caller) { 
    sapply(results.binned, function(r) {
      if(is.null(r))
        return(0)
      else {
        print("Calculating on: ")
        # print(r)
        calc(as.data.frame(r), sim.caller)
      }
    })
  })
  names(result.binned.sensitivity)=callers
  return(result.binned.sensitivity)
}

calc.binned.cnv.truth.set = function(results, truth, bin.by, bin.levels, calc=cnv.caller.recall) {
  
  if(class(results) != "GRanges")
    stop(sprintf("Results should be combined result set loaded by load.combined.results as GRanges object, but it has class %s"), class(results))
  
  truth.bins = cut(bin.by, bin.levels)
  results.binned = bin.cnv.results(truth, truth.bins, results)
  
  # Compute sensitivity for each bin
  result.binned.sensitivity = lapply(sim.callers, function(sim.caller) { 
    sapply(results.binned, function(r) {
      if(is.null(r))
        return(0)
      else {
        print("Calculating on: ")
        print(r)
        calc(as.data.frame(r), sim.caller)
      }
    })
  })
  names(result.binned.sensitivity)=sim.callers
  return(result.binned.sensitivity)
}



