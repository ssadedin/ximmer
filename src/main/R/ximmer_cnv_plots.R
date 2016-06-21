
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

XIMMER_CALLERS=unlist(strsplit(Sys.getenv("XIMMER_CALLERS"),split=","))

ANALYSIS=Sys.getenv("ANALYSIS")
if(is.na(ANALYSIS))
    ANALYSIS = "analysis"

print(sprintf("SRC=%s", SRC))

# TARGET_REGION="/Users/simon/work/ximmer/eval/data/haloplex/target_regions.bed"
TARGET_REGION=Sys.getenv("TARGET_REGION")

source(sprintf("%s/cnv_utils.R", SRC))
source(sprintf("%s/cnv_roc_lib.R", SRC))

sim.callers=XIMMER_CALLERS

sim.names = paste0("run",0:(XIMMER_RUNS-1))
ximmer.sims = lapply(sim.names, function(run) { list(name=run)})
names(ximmer.sims) = sim.names

truth = load.truth(names(ximmer.sims))
truth$id = paste(truth$id, truth$source, sep="_")

combined.results = load.combined.results(ximmer.sims, ANALYSIS, pattern=paste0("%s/",ANALYSIS,"/report/cnv_report.tsv"))
combined.results$sample = paste(combined.results$sample,combined.results$sim,sep='_')


ximmer_sim_loaders = list(
  
  ed=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) load_exomedepth_results(sprintf("%s/%s/exome_depth/%s.exome_depth.cnvs.tsv", sim, ANALYSIS, ANALYSIS),sim)))
  },
  
  xhmm=function(batch,sims) {

    do.call(c,lapply(sims, {
      function(sim) {
        glob = sprintf("%s/%s/xhmm/*.xhmm_discover.xcnv",sim, ANALYSIS)
        print(sprintf("glob = %s",glob))  
        load_xhmm_results(Sys.glob(glob)[[1]],sim)
      } 
    }))
  },
  
  cnmops=function(batch,sims) {
    do.call(c,lapply(sims,
       function(sim) {
         load_cn_mops_results(Sys.glob(sprintf("%s/%s/cn_mops/*.cnmops.cnvs.tsv",sim, ANALYSIS)[[1]]), sim)
       }))
  },
  
  cfr=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) load_conifer_results(Sys.glob(sprintf("%s/%s/conifer/*.conifer.cnvs.tsv",sim, ANALYSIS)[[1]]),sim)))
  },
  
  anghmm=function(batch,sims) {
    do.call(c,lapply(sims,function(sim) {
      pattern = sprintf("%s/%s/%s.*counts.angelhmm*.cnvs.bed", sim, ANALYSIS, ANALYSIS)
      print(sprintf("Pattern = %s", pattern))
      load_angel_results(Sys.glob(pattern)[[1]],sim)
    }))
  }
)

#ximmer.sims = list(
#  "run0" = list(name = "run0"),
#  "run1" = list(name = "run1"),
#  "run2" = list(name = "run2")
#)

sim.caller.loaders = ximmer_sim_loaders
ranked = load_ranked_run_results(truth, names(ximmer.sims), filterChrX=F)
nice_colors = plot.colors = c('orange','blue','green','black','purple','red')

png("roc.png")
plot.exome.result.set(ranked, truth)
dev.off()

truth.metrics = compute_deletion_metrics(truth, TARGET_REGION)

png("sens_by_tg.png")
bin.levels = c(0,1,2,3,4,5,6,7,8)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$targets, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=8, plot.xlab = "Deletion Size (No. of Target Regions)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png("sens_by_seqbp.png")
bin.levels = c(0,200,500,1000,1500,4000)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$seqbp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=4000, plot.xlab = "Deletion Size (sequenced bp)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png("sens_by_bp.png")
bin.levels = c(0,200,500,1000,1500,4000,10000,100000)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$bp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, plot.xlab = "Deletion Size (spanned bp)", log.scale = T)
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors, log.scale = T)
dev.off()

##--------------- Quality Score Calibration

png("qual_score_calibration.png", width=960, height=250*length(sim.callers))
par(mfrow=c(length(sim.callers)/2,2))
for(caller in sim.callers) {
  if(length(ranked[[caller]])> 0)
    plot.qscore.calibration(caller, ranked, newPlot=T, col="darkgreen")
}
dev.off()


