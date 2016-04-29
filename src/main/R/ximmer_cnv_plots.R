
# setwd("/Users/simon/work/ximmer/debug_runs")

TOOLS = Sys.getenv("TOOLS")
if(is.na(TOOLS) || TOOLS=="")
  TOOLS="/Users/simon/work/ximmer/eval/tools"

SRC = Sys.getenv("SRC")
if(is.na(SRC) || SRC=="")
  SRC="/Users/simon/work/ximmer/src/main/R"

XIMMER_RUNS = as.integer(Sys.getenv("XIMMER_RUNS"))
if(is.na(XIMMER_RUNS))
    XIMMER_RUNS = 5

# TARGET_REGION="/Users/simon/work/ximmer/eval/data/haloplex/target_regions.bed"
TARGET_REGION=Sys.getenv("TARGET_REGION")

print(sprintf("SRC=%s", SRC))

source(sprintf("%s/cnv_utils.R", SRC))
source(sprintf("%s/cnv_roc_lib.R", SRC))

ximmer_sim_loaders = list(
  
  ed=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) load_exomedepth_results(sprintf("%s/analysis/analysis.exome_depth.cnvs.tsv", sim),sim)))
  },
  
  xhmm=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) load_xhmm_results(Sys.glob(sprintf("%s/analysis/xhmm/*.xhmm_discover.xcnv",sim)[[1]]),sim)))
  },
  
  mops=function(batch,sims) {
    do.call(c,lapply(sims,
       function(sim) {
         load_cn_mops_results(Sys.glob(sprintf("%s/analysis/cn_mops/*.cnmops.cnvs.tsv",sim)[[1]]), sim)
       }))
  },
  
  cfr=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) load_conifer_results(Sys.glob(sprintf("%s/analysis/*.conifer.cnvs.tsv",sim)[[1]]),sim)))
  },
  
  anghmm=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) {
      pattern = sprintf("%s/analysis/analysis.*counts.angelhmm*.cnvs.bed", sim)
      print(sprintf("Pattern = %s", pattern))
      load_angel_results(Sys.glob(pattern)[[1]],sim)
    }))
  }
)

ximmer.sims = list(
  "run0" = list(name = "run0"),
  "run1" = list(name = "run1"),
  "run2" = list(name = "run2")
)

sim.names = paste0("run",0:(XIMMER_RUNS-1))
ximmer.sims = lapply(sim.names, function(run) { list(name=run)})
names(ximmer.sims) = sim.names

truth = load.truth(names(ximmer.sims))
truth$id = paste(truth$id, truth$source, sep="_")

sim.callers=c("mops","xhmm")
sim.caller.loaders = ximmer_sim_loaders
ranked = load_ranked_run_results(truth, names(ximmer.sims), filterChrX=F)
nice_colors = plot.colors = c('orange','blue','green','black','purple','red')

png("roc.png")
plot.exome.result.set(ranked, truth)
dev.off()

truth.metrics = compute_deletion_metrics(truth, TARGET_REGION)
x = load.combined.results(ximmer.sims, "analysis", pattern="%s/analysis/report/cnv_report.tsv")
x$sample = paste(x$sample,x$sim,sep='_')

png("sens_by_tg.png")
bin.levels = c(0,1,2,3,4,5,6,7,8)
x.binned = load.binned.cnv.truth.set(x, "analysis", truth.metrics, truth.metrics$targets, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=8, plot.xlab = "Deletion Size (No. of Target Regions)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png("sens_by_seqbp.png")
bin.levels = c(0,200,500,1000,1500,4000)
x.binned = load.binned.cnv.truth.set(x, "analysis", truth.metrics, truth.metrics$seqbp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=4000, plot.xlab = "Deletion Size (sequenced bp)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png("sens_by_bp.png")
bin.levels = c(0,200,500,1000,1500,4000,10000,100000)
x.binned = load.binned.cnv.truth.set(x, "analysis", truth.metrics, truth.metrics$bp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, plot.xlab = "Deletion Size (spanned bp)", log.scale = T)
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors, log.scale = T)
dev.off()



 
