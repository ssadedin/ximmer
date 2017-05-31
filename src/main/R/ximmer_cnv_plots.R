
#setwd("/Users/simon/work/ximmer/debug_runs")

TOOLS = Sys.getenv("TOOLS")
if(is.na(TOOLS) || TOOLS=="")
  TOOLS="/Users/simon/work/ximmer/eval/tools"

SRC = Sys.getenv("SRC")
if(is.na(SRC) || SRC=="")
  SRC="/Users/simon/work/ximmer/src/main/R"

XIMMER_RUNS = unlist(strsplit(Sys.getenv("XIMMER_RUNS"), ","))
if(is.na(XIMMER_RUNS))
    XIMMER_RUNS = c("1")

# XIMMER_CALLERS=c('cnmops_1e1','cnmops_1e2')
XIMMER_CALLERS=unlist(strsplit(Sys.getenv("XIMMER_CALLERS"),split=","))

ANALYSIS=Sys.getenv("ANALYSIS")
if(is.na(ANALYSIS))
    ANALYSIS = "analysis-cnmops_tuning"

SIMULATION_TYPE=Sys.getenv("SIMULATION_TYPE")
if(is.na(SIMULATION_TYPE))
    SIMULATION_TYPE = "replace"


print(sprintf("SRC=%s", SRC))
print(sprintf("Callers=%s", paste(XIMMER_CALLERS,collapse=',')))

#TARGET_REGION="/Users/simon/work/ximmer/eval/data/haloplex/target_regions.bed"
#TARGET_REGION="/Users/simon/work/ximmer/eval/data/haloplex.x/target_regions.bed"
TARGET_REGION=Sys.getenv("TARGET_REGION")

source(sprintf("%s/cnv_utils.R", SRC))
source(sprintf("%s/cnv_roc_lib.R", SRC))

sim.callers=XIMMER_CALLERS

# sim.callers = c("xhmm_1","xhmm_2")

sim.names = XIMMER_RUNS
ximmer.sims = lapply(sim.names, function(run) { 
  list(name=run)}
)

# XIMMER_CALLER_LABELS = c('pi=5','pi=10')
XIMMER_CALLER_LABELS = unlist(strsplit(Sys.getenv('XIMMER_CALLER_LABELS'), split=','))

print(sprintf("Labels = %s", paste(XIMMER_CALLER_LABELS,sep=',',collapse=',')))

# Load DGV

DGV_MAX_FREQ=as.numeric(Sys.getenv("DGV_MAX_FREQ"))
DGV_MIN_STUDY_SIZE=as.numeric(Sys.getenv("DGV_MIN_STUDY_SIZE"))
DGV_CNVS=Sys.getenv("DGV_CNVS")

print(sprintf("Loading cnvs with min freq %f, min study size %s from %s", 
            DGV_MAX_FREQ, DGV_MIN_STUDY_SIZE, DGV_CNVS))
       
dgv = load_dgv(DGV_CNVS)

print(sprintf("Filtering DGV CNVs"))

dgv.exclude = dgv[(dgv$count > DGV_MIN_STUDY_SIZE) & (dgv$perc>DGV_MAX_FREQ)]

# print(sprintf("Excluding %d regions from DGV", sum(width(reduce(dgv.exclude)))));

print("Loaded DGV CNVs");

names(ximmer.sims) = sim.names

truth = load.truth(names(ximmer.sims))
truth$id = paste(truth$id, truth$source, sep="_")

combined.results = load.combined.results(ximmer.sims, ANALYSIS, pattern=paste0("%s/",ANALYSIS,"/report/cnv_report.tsv"))
combined.results$sample = paste(combined.results$sample,combined.results$sim,sep='_')

combined.results

ximmer_sim_loaders = list(
  
  ed=function(label, batch,sims) {
    do.call(c,lapply(sims,function(sim) load_exomedepth_results(sprintf("%s/%s/%s/%s.exome_depth.cnvs.tsv", sim, ANALYSIS, label, ANALYSIS),sim)))
  },
  
  xhmm=function(label, batch,sims) {

    do.call(c,lapply(sims, {
      function(sim) {
        glob = sprintf("%s/%s/%s/*.xhmm_discover.xcnv",sim, ANALYSIS,label)
        print(sprintf("glob = %s",glob))  
        load_xhmm_results(Sys.glob(glob)[[1]],sim)
      } 
    }))
  },
  
  cnmops=function(label, batch,sims) {
    do.call(c,lapply(sims,
       function(sim) {
         load_cn_mops_results(Sys.glob(sprintf("%s/%s/%s/*.cnmops.cnvs.tsv",sim, ANALYSIS,label)[[1]]), sim)
       }))
  },
  
  cfr=function(label, batch,sims) {
    do.call(c,lapply(sims,function(sim) load_conifer_results(Sys.glob(sprintf("%s/%s/%s/*.conifer.cnvs.tsv",sim, ANALYSIS, label)[[1]]),sim)))
  },
  
  anghmm=function(label, batch,sims) {
    do.call(c,lapply(sims,function(sim) {
      pattern = sprintf("%s/%s/%s/%s.*counts.angelhmm*.cnvs.bed", sim, ANALYSIS, label, ANALYSIS) # note: not tested yet for angel, may need pipeline adjustment
      print(sprintf("Pattern = %s", pattern))
      load_angel_results(Sys.glob(pattern)[[1]],sim)
    }))
  }
)

sim.caller.loaders = ximmer_sim_loaders

filterChrX=FALSE
deletionsOnly=FALSE

if(SIMULATION_TYPE == "replace") {
  filterChrX=TRUE
  deletionsOnly=TRUE
} else if(SIMULATION_TYPE == "downsample") {
    deletionsOnly=TRUE
}

ranked = load_ranked_run_results(truth, names(ximmer.sims), exclude.bed=dgv.exclude, filterChrX=filterChrX, deletionsOnly=deletionsOnly)
nice_colors = plot.colors = c('orange','blue','green','black','purple','red', 'pink', 'brown','gray','aquamarine','chartreuse','coral')

sim.caller.labels = XIMMER_CALLER_LABELS
names(sim.caller.labels) = XIMMER_CALLERS

png(sprintf("%s_roc.png",ANALYSIS))
plot.exome.result.set(ranked, truth)
dev.off()

truth.metrics = compute_deletion_metrics(truth, TARGET_REGION)

png(sprintf("%s_sens_by_tg.png", ANALYSIS))
bin.levels = c(0,1,2,3,4,5,6,7,8)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$targets, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=8, plot.xlab = "Deletion Size (No. of Target Regions)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png(sprintf("%s_sens_by_seqbp.png", ANALYSIS))
bin.levels = c(0,200,500,1000,1500,4000)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$seqbp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, x.max=4000, plot.xlab = "Deletion Size (sequenced bp)")
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors)
dev.off()

png(sprintf("%s_sens_by_bp.png", ANALYSIS))
bin.levels = c(0,200,500,1000,1500,4000,10000,100000)
x.binned = load.binned.cnv.truth.set(combined.results, ANALYSIS, truth.metrics, truth.metrics$bp, bin.levels)
sens_by_del_plot_frame(bin.levels = bin.levels, plot.xlab = "Deletion Size (spanned bp)", log.scale = T)
plot.binned.performance(r=x.binned, bin.levels=bin.levels, palette=plot.colors, log.scale = T)
dev.off()

##--------------- Quality Score Calibration

png(sprintf("%s_qual_score_calibration.png",ANALYSIS), width=960, height=250*length(sim.callers))
if(length(sim.callers)>1) {
  par(mfrow=c(length(sim.callers)/2,2))
}
for(caller in sim.callers) {
  if(length(ranked[[caller]])> 0)
    plot.qscore.calibration(caller, ranked, newPlot=T, col="darkgreen")
}
dev.off()


