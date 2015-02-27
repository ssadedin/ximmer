library(VariantAnnotation)
library("rootSolve")
library(MASS)

hg19.chromosomes = paste0("chr",c(1:22, "X","Y"))

read.bed = function(f) {
  read.table(f, col.names=c("chr","start","end","id"), fill=1)
}

read.bed.ranges = function(f) {
  bed = read.table(f, col.names=c("chr","start","end","id"), fill=1)
  GRanges(seqnames=bed$chr, ranges=IRanges(start=bed$start, end=bed$end-1), id=bed$id)
}

load_xhmm_results = function(file.name,sample.tag=NA) {
  # Load output from XHMM and return it in a "standardized" format
  # CHR START END TYPE QUALITY SAMPLE
  raw = read.table(file.name, header=T, colClasses=c("character"))
  raw$start = as.integer(gsub("^.*:","",gsub("-.*$","",raw$INTERVAL)))
  raw$end = as.integer(gsub("^.*-","",raw$INTERVAL))
  
  # Trim any trailing underscores from sample name (a hack to fix some that got left there by accident once)
  if(!is.na(sample.tag) && nrow(raw))
    xhmm.samples = paste(xhmm.samples,sample.tag,sep="_")

  # Remove any entries not in hg19.chromosomes as this will create issues later
  raw = raw[raw$CHR %in% hg19.chromosomes,]
  
  xhmm.samples = gsub("_$","",raw$SAMPLE)

  return(GRanges(seqnames=raw$CHR, 
                 seqinfo=Seqinfo(hg19.chromosomes),
                 ranges=IRanges(start=raw$start, end=raw$end), type=as.character(raw$CNV), sample=xhmm.samples, qual=raw$Q_SOME))
}

load_excavator_results = function(file.name, sample.tag=NA) {
  cnv.types = c("DEL", "DEL","NA", "DUP","DUP") # Will access this by index
  raw = read.table(file.name, col.names=c("sample","chr","start","end","ratio","cn.change"),header=F, as.is=T)
  raw$sample= gsub('_$','',as.character(raw$sample))
  if(!is.na(sample.tag) && nrow(raw))
    raw$sample = paste(raw$sample,sample.tag,sep="_")
  
  # Remove any entries not in hg19.chromosomes as this will create issues later
  raw = raw[raw$chr %in% hg19.chromosomes,]
  
  result = (GRanges(seqnames=raw$chr,
                 seqinfo=Seqinfo(hg19.chromosomes),
                 ranges=IRanges(start=raw$start, end=raw$end),
    type=as.character(unlist(lapply(raw$cn.change+3, function(x) {cnv.types[[x]]}))), sample=raw$sample))
  
  if(length(result)>0)
    result$qual=log10(raw$end - raw$start)
    
  return(result);
}

load_exomedepth_results = function(file.name, sample.tag=NA) {
  cnv.types = list(deletion="DEL", duplication="DUP")
  raw = read.table(file.name, header=T, as.is=T, colClasses=c("integer","integer","character","integer","integer","integer","character","character","double","integer","integer","double","character"))
  # Some samples have trailing underscores incorrectly
  raw$sample= gsub('_$','',as.character(raw$sample))
  if(!is.na(sample.tag) && nrow(raw))
    raw$sample = paste(raw$sample,sample.tag,sep="_")
  
  raw = raw[raw$chromosome %in% hg19.chromosomes,]

  return(GRanges(seqnames=raw$chromosome, 
                 seqinfo=Seqinfo(hg19.chromosomes),
                 ranges=IRanges(start=raw$start, end=raw$end), 
                 type=unlist(lapply(raw$type, function(x) {cnv.types[[x]]})), sample=raw$sample, qual=raw$BF))
}

load_angel_results = function(file.name, sample.tag=NA) { 
  # Load results from the Angel output format, which is BED-like but
  # has an extra column giving the quality estimate for each deletion call
  
  raw = read.table(file.name, 
             as.is=T,
             col.names=c("chr","start","end","sample","qual"),fill=1) # quality not present for truth file
  
  if(length(raw)>0)
	  raw$type = 'DEL'
  
  # NOTE: replace trailing _ because of some earlier mistakes
  # with them not being trimmed off and ending up in the BAM headers
  raw$sample= gsub('_$','',as.character(raw$sample))
  
  if(!is.na(sample.tag))
      raw$sample = paste(raw$sample,sample.tag,sep="_")
  
  raw = raw[raw$chr %in% hg19.chromosomes,]

  return(GRanges(seqnames=raw$chr,
                 seqinfo=Seqinfo(hg19.chromosomes),
                 ranges=IRanges(start=raw$start, end=raw$end),
                 type=raw$type, sample=raw$sample, qual=raw$qual))
}

# The "truth" results are saved in the same format as the angel results
load_truth_results = load_angel_results

load_exome_copy_results = function(file.name, sample.tag=NA) {
  raw = read.table(file.name, header = T)  
  
  if(!is.na(sample.tag))
      raw$sample.name = paste(raw$sample.name,sample.tag,sep="_")
    
  raw = raw[raw$space %in% hg19.chromosomes,]

  result = GRanges(seqnames=raw$space, 
                   ranges=IRanges(start=raw$start, end=raw$end),
                   type=ifelse(raw$copy.count<2,'DEL','DUP'),
                   sample=raw$sample.name,
                   qual=raw$log.odds
                   )
  return(result)
}

load_cn_mops_results = function(file.name, sample.tag=NA) {
  raw = read.table(file.name, header = T, colClasses=c(NA,NA,NA,NA,NA,"character"))  
  raw$CN =  as.integer(gsub("^CN","",raw$CN))
  raw = raw[raw$CN != 2,]
  
  if(!is.na(sample.tag))
      raw$sampleName = paste(raw$sampleName,sample.tag,sep="_")
    
  raw = raw[raw$seqnames %in% hg19.chromosomes,]

  result = GRanges(seqnames=raw$seqnames, 
                   ranges=IRanges(start=raw$start, end=raw$end),
                   type=ifelse(raw$CN<2,'DEL','DUP'),
                   sample=raw$sampleName,
                   qual=-raw$median
                   )
  return(result)
}


load_dgv = function(gzipped.file.name) {
  # Read the UCSC database of genomic variants file
  # into a GRanges object and set a few common fields on there
  dgv = read.table(pipe(sprintf("gunzip -c %s", gzipped.file.name)), header=F, fill=5, sep='\t')
  names(dgv)=c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "varType", "reference", "pubMedId", "method", "platform", "mergedVariants", "supportingVariants", "sampleSize", "observedGains", "observedLosses", "cohortDescription", "genes", "samples")
  
  dgv.ranges = GRanges(seqnames=dgv$chrom, ranges=IRanges(start=dgv$chromStart, end=dgv$chromEnd))
  dgv.ranges$count = dgv$sampleSize
  dgv.ranges$gains = dgv$observedGains
  dgv.ranges$losses = dgv$observedLosses
  return(dgv)
}


#####################################################################################

combine_cnv_caller_results = function(cnv.samples, cnv.results, target.bed, chr="*") {
  #
  # Combine separate sets of calls from the different CNV callers
  # and return a list indexed by sample id, containing the aggregated
  # results from all the callers for each sample separately
  #
  # Arguments:
  #
  #    cnv.samples - a list of samples to aggregate the results for
  #    cnv.results - a list of results for different CNV callers,
  #                  indexed by the name of the caller (ex, ed, xhmm, etc.)
  #    target.bed  - a GRanges containing the target region with
  #                  gene annotations loaded as the 'id' column
  #                  (see read.bed.ranges function)
  
  
  # Compute merged set on a per sample basis
  cnv.all=GRanges()
  cnv.all.merge = list()
  seqlevels(cnv.all)=hg19.chromosomes
  for(s in cnv.samples) {
    print(s)
    merged = extract_merged_cnvs(s, cnv.results, chr)
    cnv.all.merge[[s]] = merged
    cnv.all = c(cnv.all, merged)  
  }
  
  #cnv.all$count = rowSums(as.matrix(mcols(cnv.all)[,c('xhmm','ed','ex')]))
  cnv.all$count = rowSums(as.matrix(mcols(cnv.all)[,names(cnv.results)]))
  
  # compute a "hash" of each unique CNV to make it easy to plot a venn diagram
  cnv.all$chrpos=paste(cnv.all$sample,paste(paste(seqnames(cnv.all),start(cnv.all),sep=":"),end(cnv.all),sep="-"),sep="/")
  
  # add gene annotations to the cnvs from the target bed file
  cnv.all.target.hits = as.data.frame(findOverlaps(cnv.all,target.bed))  
  
  if(length(cnv.all)>0) {
    cnv.all$genes = aggregate(cnv.all.target.hits$subjectHits, list(cnv.all.target.hits$queryHits), 
                        function(target.hits) { paste(unique(target.bed[target.hits]$id),collapse=',')})$x
  }
  
  return(cnv.all)
}

extract_merged_cnvs = function(s, cnv.results, chr="*") {
  #
  # Extracts the CNVs for a particular sample from the list of CNV results given
  # that are called by separate callers. The results are expected to be in the 
  # standard format as loaded by the load_xxx_results functions above.
  #
  # Args:
  #    -  cnv.results - a list of GRanges objects, one for each caller
  #
  #    -            s - the sample of interest
  #
  
  # Extract the CNVs belonging to sample s from each caller
  s.cnvs = sapply(cnv.results, function(cnv.caller.results) {
    if(chr=="*")
      cnv.caller.results[cnv.caller.results$sample==s]
    else
      cnv.caller.results[cnv.caller.results$sample==s & (seqnames(cnv.caller.results)==chr)]
  },simplify=F, USE.NAMES=T)
  
  merged = GRanges()
  seqlevels(merged)=hg19.chromosomes;
  for(xx in s.cnvs) { seqlevels(xx) = hg19.chromosomes;merged = c(merged,xx) }
  
  merged = reduce(merged)

  for(caller in names(s.cnvs)) {
    mcols(merged)[,caller]=merged %over% s.cnvs[[caller]]
  }
  
  if(length(merged)>0) {
    merged$sample = s
    merged$type="" # to be replaced below
    for(i in 1:length(merged)) {
      types = unlist(lapply(s.cnvs, function(cnvs) { 
        s.overlaps = cnvs[subjectHits(findOverlaps(merged[i], cnvs))]
        if(length(s.overlaps)) 
          as.vector(s.overlaps$type)[[1]]
        else
          c()
      }))
      if(length(unique(types)) != 1) {
        print("Warning: cnv region has conflicting calls")
        print(merged[i])
        print(types)
      }
      merged[i]$type = types[[1]]
    }
  }  
  return(merged);
}

#####################################################################################

find.cnv.overlaps = function(region, cnvs, sample=F) {
  if(sample == F) {
    sample=region$sample[[1]]
  }
  cnvs[cnvs$sample==sample][subjectHits(findOverlaps(region,cnvs[cnvs$sample==sample]))]
}

#####################################################################################

load.normalised.coverage = function(samples, cov.files, scale="relative", chr="*") {
  # Load coverage values as computed by coverageBed and normalise to 
  # account for number of reads per sample, and then 
  # per position to give values relative to mean coverage for all samples
  # at each position.
  #
  # Args:
  #   samples   list of sample names to load coverage for
  #
  #   cov.files either a list with names corresponding to the sample names containing the
  #             coverage file for each sample OR a pattern names containing %s into which 
  #             to substitute the sample name to obtain the file name for each sample
  #  
  #   scale     either "relative" or "absolute". If relative, the plot
  #             will be scaled so that values are first scaled relative to their means
  #             and then relative to the mean of all samples at each position so that the
  #             "expected" coverage level is at 1.0.
  #             If "absolute", the values will be normalised relative to their 
  #             means, but then scaled up to be relative to the median of all samples.
  #             The purpose of this is to give a feeling for the absolute coverage at
  #             each position.
  #          
  
  if(typeof(cov.files)=="list") {
      print("Coverage files provided explicitly")
      names(cov.files)=samples
  }
  else {
      cov.files = sapply(samples, function(s) {  Sys.glob(sprintf(cov.files,s)) }, USE.NAMES=T)
  }
  
  awk_filter = ""
  if(chr!="*") {
    awk_filter = sprintf(" | awk '{ if($1 == \"%s\") print $0 }' ", chr)
  }
  cmd = sprintf(paste0("cut -f 1,2,3,4,5 ",cov.files[[1]], awk_filter),samples[[1]])
  print(cmd)
  cov.info = read.table(pipe(cmd),
                            col.names=c("chr","exon.start","exon.end","gene","exon.pos"), header=F)
  
  cov.info$pos = cov.info$exon.start+cov.info$exon.pos-1
  
  print(sprintf("Number of rows = %d", nrow(cov.info)))
  
  cov.tmp = matrix(nrow=nrow(cov.info))
  for(s in samples) {
    print(s)
    #cov.tmp = cbind(cov.tmp, read.table(pipe(fmt("cut -f 6 <%=base.dir%>/work/cnv/<%=s%>_<%=machine.name%>_*.exoncoverage.txt")),header=F)$V1)
    cmd=sprintf(paste0("cut -f 1,6 ",cov.files[[s]],awk_filter," | cut -f 2 "),s)
    print(sprintf("Command = %s", cmd))
    cov.tmp = cbind(cov.tmp, read.table(pipe(cmd),header=F)$V1)
  }
  cov.tmp= cov.tmp[,-1, drop=F] # Note first column is NA's, chop it off
  colnames(cov.tmp) = samples
  
  # Normalise by dividing out the mean coverage for each sample
  normcov = cov.tmp / colMeans(cov.tmp)

  cov.scaled = normcov * mean(rowMeans(cov.tmp))
  
  # Normalise by position by dividing out mean coverage of each position
  # so each sample is now relative to the other samples at that position
  
  #browser()
  normcov = normcov / rowMeans(normcov)
  
  normcov[is.nan(normcov)] <- 0
  
  #normcov = normcov / ifelse(means>0,means,1)
  
  return(list(info=cov.info, cov.abs=cov.tmp, cov.norm=normcov, cov.scaled=cov.scaled ))
}

#####################################################################################


plot.hbar = function(x.pos, y.pos, label, col='blue',label.pos='above', cex.text=0.7, height=0.03, ...) {
  lines(x.pos,rep(y.pos,2), col=col, ...) 
  lines(c(x.pos[[1]],x.pos[[1]]), c(y.pos-height,y.pos+height),col=col, ...)
  lines(c(x.pos[[2]],x.pos[[2]]), c(y.pos-height,y.pos+height),col=col, ...)
  if(label.pos=='above')
    text(mean(x.pos), y.pos+strheight(label,cex=cex.text)*1.4, label,col=col,cex=cex.text)
  else
  if(label.pos=='left')
    text(x.pos[[1]]-strwidth(label,cex=cex.text), y.pos, label,col=col,cex=cex.text)
  else
    stop(paste0("Invalid valid of label.pos parameter: ",label.pos))
    
}

plot.cnv = function(cnv, all.cnvs, cov.data, variants=NA, add.title.text=T, cex.cnv.label=0.7, scale="relative",...) {
  # 
  # Finds all CNVs that overlap the region of a particular CNV and then 
  # plots them all together on the same plot.
  #
  # Args: 
  #   cnv       : a GRanges object including 'sample' column indicating the 
  #               CNV (or range) and sample to use for extracting and plotting 
  #               CNVs. Other CNVs will be plotted too if they appear in the 
  #               same region.
  #
  #   all.cnvs :  A list of GRanges objects containing all CNVs to be considered. 
  #               Each GRanges object represents CNV calls from a single CNV caller. 
  #               The names of the list are assumed to be names of CNV callers and will
  #               be used to label the CNVs according to which CNV caller 
  #               called them.
  #  
  #   variants :  a data frame containing columns pos, gt, dp, alt.dp, ref.dp
  #               see find.cnv.variants for how to create this.
  #
  #
  results = list()
  for(caller in names(all.cnvs)) {
    results[[caller]] = find.cnv.overlaps(cnv, all.cnvs[[caller]])    
  }
  if(typeof(variants)=="S4" | typeof(variants)=="list") {
    if(typeof(variants)=="list") {
      variants = variants[[cnv$sample]]  
    }
    variants_to_plot = find.cnv.variants(variants, cnv)
  }
  else {
    variants_to_plot = FALSE
  }
  plot.cnv.region(results, cov.data, variants_to_plot, add.title.text=add.title.text, cex.cnv.label=cex.cnv.label, scale=scale,...)
}

plot.cnv.region = function(cnvs, 
                          cov.data, 
                          variants=F, 
                          highlight.amplicons=F, 
                          add.title.text=T, 
                          cex.cnv.label=0.7, 
                          xlim=NULL,
                          scale="relative",
                          show.all.samples=F,...) {
  #
  # Plots coverage over a region containing a CNV, showing mean coverage
  # from all samples vs a single sample containing a CNV
  #
  # Args:
  #    cnvs     : GRanges object defining region to plot, or list of them with names
  #               The ranges must have a meta data value called 'sample' indicating the 
  #               sample name for the CNV.
  #  
  #    cov.data : Coverage data for the targeted sequencing region, as loaded by 
  #               "load_normalised_coverage" function, ie: containing 
  #               attributes "info" and "cov.norm".
  #
  #    variants : A data frame of variants may be provided, which are then plotted
  #               at their respective positions with indicator of het/hom
  #               and number of reads
  #
  #    scale    : If "relative", coverage will be scaled such that it is relative to
  #               the mean of all samples at each base position, thus making 1.0 the
  #               straight line representing "expected" coverage and the sample of interest
  #               as a fraction of that coverage level. If "absolute", 
  #               the expected coverage will be represented as the median value of normalised
  #               coverage (where the normalisation is done by dividing out the mean from each sample,
  #               but then scaling coverage back up to the mean of all samples.)
  #               
  #
  cov.norm = cov.data$cov.norm
  cov.info = cov.data$info
  cov.abs = cov.data$cov.abs[order(cov.info$chr,cov.info$pos),]
  
  if(typeof(cnvs)=="list") {
    cnv.groups=cnvs
  } else {
    cnv.groups=list(cnv=cnvs)  
  }
  
  # Extract sample name from the CNVs to plot
  s = unique(unlist(lapply(cnv.groups, function(cr) { unique(cr$sample)})))
  if(length(s)>1)
    stop("List of CNVs given contains CNVs from more than one sample")
  
  other.samples = as.character(colnames(cov.norm)[colnames(cov.norm) != s])
  print(sprintf("Other samples are %s", paste(other.samples,collapse=",")))
  
  cov.scaled = (cov.abs / colMeans(cov.abs)) * mean(colMeans(cov.abs))

  smoothing = 0.003

  cnvs.tmp = cnv.groups
  names(cnvs.tmp)=NULL # remove names
  #browser()
  region = reduce(do.call("c", cnvs.tmp))
  
  #print(fmt("Region is <%=as.character(seqnames(region))%>:<%=start(region)%>-<%=end(region)%>"))
  
  # Sort the coverage information

  cov.norm = cov.norm[order(cov.info$chr,cov.info$pos),]
  cov.info = cov.info[order(cov.info$chr,cov.info$pos),]
  
  # Find all the genes covered by the region
  region.genes = cov.info[cov.info$chr==as.character(seqnames(region)) & cov.info$pos>start(region) & cov.info$pos<end(region),'gene']
  
  start.gene = region.genes[[1]]
  end.gene = region.genes[[length(region.genes)]]
  
  # Find the region from start of the first gene, to end of last gene
  region.info = cov.info[min(match(start.gene,cov.info$gene)):(nrow(cov.info)-match(end.gene,rev(cov.info$gene))),]
  
  # Show no more than 1kb downstream and upstream of the start and end gene
  plot.start = max(0,which(region.info$pos == start(region))-1000 )
  plot.end = min(which(region.info$pos == end(region))+1000, nrow(region.info))
  
  region.info = region.info[plot.start:plot.end,]
  
  region.positions = cov.info$chr==region.info[1,'chr'] & cov.info$pos>region.info[1,'pos'] & cov.info$pos<region.info[nrow(region.info),'pos']
  
  region.scaled.cov = cov.scaled[region.positions,]
  
  region.other.scaled.meds = apply(region.scaled.cov[,other.samples], 1, mean)
  region.other.scaled.meds.log = log2(apply(region.scaled.cov[,other.samples], 1, mean))
    
  region.cov = cov.norm[region.positions,]
  
    # Draw the general coverage over the gene
  par(xpd=F)
  par(mar=c(12,5,5,5))
  if(scale == "relative") {
    y.max = 2.0
    cov.plot = region.cov
    cov.plot.others = region.cov
    ylab="Relative Normalised Coverage"
  }
  else {
    y.max = max(lowess(region.other.scaled.meds.log, f=smoothing)$y) * 1.10
    #browser()
    cov.plot = log2(region.scaled.cov)
    ylab="Log of Normalised Absolute Coverage"
  }
  print(sprintf("Y max = %f",y.max))
  
  
#  plot(lowess(cov.plot[,s],f=smoothing), t="l", ylim=c(0,y.max), 
#       #xlab="Position in Region", 
#       xlab="",
#       ylab=ylab,
#       main="",...) #, xlim=c(1000,3000))

  pos_to_x = function(pos) {
    which.min(abs(pos-region.info$pos))
  }

  exon.starts = unique(region.info$exon.start)
  exon.ends = unique(region.info$exon.end)
  
  exon.starts.x = sapply(exon.starts, pos_to_x)
  exon.ends.x = sapply(exon.ends, pos_to_x)

  if(missing(xlim)) {
    xlim = c(0,nrow(region.info))
  }
  
  plot(c(), t="n", ylim=c(0,y.max), 
       #plot(cov.plot[,s], t="l", ylim=c(0,y.max), 
       #xlab="Position in Region", 
       xlab="",
       ylab="Relative Normalised Coverage",
       main="", xlim=xlim,...)
  

  y_per_inch = y.max / par()$pin[2]
  x_per_inch = nrow(cov.plot) / par()$pin[1]
  #print(paste0("Y per inch = ",y_per_inch))
  
  title(xlab="Position in Region", line=1.9)
  
  # extra smoothed line?
  #lines(lowess(region.cov[,s],f=smoothing*10), t="l", lty=2, ylim=c(0,2.0), col="#555555aa")
    
  if(scale=="relative") {
    #lines(lowess(1+apply(region.cov,1,sd), f=smoothing), col=2, lty=2)
    sdlowess=lowess(1-apply(region.cov[,other.samples],1,sd), f=smoothing)
    lines(sdlowess, col=2, lty=2)
    polygon(c(1:nrow(region.info), rev(sdlowess$x)),
            c(rep(1,nrow(region.info)), rev(sdlowess$y)),
            col="#f0e0e050",border=F
            )
    
    # Line at 0.5 to show expected coverage for het deletion
    abline(h=0.5, lty=2, col='gray')
  }
  else {
    sd.diff = (region.other.scaled.meds - apply(region.scaled.cov[,other.samples], 1, function(x) sd(x)))+1
    sd.diff.log = sapply(sd.diff, function(x) if(x>0.1) log2(x) else 0)
    lines(lowess(sd.diff.log, f=smoothing), col=2, lty=2) 
  }

  if(scale=="relative") {
    abline(h=1.0, col=2)
  }
  else {
    lines(lowess(region.other.scaled.meds.log, f=smoothing), col=2)
    #lines(region.other.scaled.meds, col=2, lty=2)
  }

  for(i in 1:length(exon.starts)) {
    starts = exon.starts.x[[i]]
    ends = min(exon.ends.x[[i]], nrow(cov.plot))
    if((length(ends)>0) && (length(starts)>0) && (ends != starts)) {
      #if(i == 3)
      #    browser()
      if(ends-starts>1) {
        if(show.all.samples) {
          for(j in other.samples) {
              lines(lowess(starts:ends,cov.plot[starts:ends,j], f=(75 * 4000 / (xlim[2]-xlim[1]))/ (ends-starts)), lwd=0.5,col='#aaaaaa',...)
          }
        }
        lines(lowess(starts:ends,cov.plot[starts:ends,s], f=(75 * 4000 / (xlim[2]-xlim[1]))/ (ends-starts)), lwd=2,...)      
      }
    }
  }
  #points(1-apply(region.cov,1,sd), col=2, lty=2, pch=19)
  
  # Draw yellow rectangle over region of actual call
  #rect(match(start(region),region.info$pos),-1,match(end(region),region.info$pos), 2.5, col="#eeee9944")
  
  text_height_inches = par()$ps / 72
    
  hbar.height = y_per_inch * 0.1
  
  # Draws the CNV calls for each caller over the region
  for(j in 1:length(cnv.groups)) {
    cnvs = cnv.groups[[j]] 
    group.name=names(cnv.groups)[[j]]
    par(xpd=NA)
    cnvs = cnv.groups[[j]]
    if(length(cnvs)>0) {
      for(i in 1:length(cnvs)) {
        print(cnvs[i])
        print(sprintf("Caller: %s (%d/%d), group %d", group.name, i, length(cnvs),j))
        x.pos = c(pos_to_x(start(cnvs[i])),pos_to_x(end(cnvs[i])))
        y.pos=y.max+y_per_inch*(0.3 + (j-1)/4) # 1/3 inch per cnv caller
        plot.hbar(x.pos, y.pos, group.name, j+1, cex.text=cex.cnv.label, label.pos='left', height=hbar.height)
      }
    }
  }
  chr=region.info$chr[[1]]
  len=length(unique(region.genes))
  
  title_text=sprintf("Sample %s Chromosome %s (%d Genes)",s,chr,len)
  title_width=strwidth(title_text,cex=1.2)
  if(add.title.text) {
      text((nrow(region.info)-title_width)/2, y.pos+y_per_inch*(0.35+text_height_inches),
         title_text,
         adj=0,cex=1.5)
  }
    
  # Draw the genes
  gene.starts = aggregate(region.info$pos, list(region.info$gene), min)
  gene.ends = aggregate(region.info$pos, list(region.info$gene), max)
  for(i in 1:nrow(gene.starts)) {
    plot.hbar(c(match(gene.starts[i,'x'],region.info$pos),
                match(gene.ends[i,'x'],region.info$pos)), 
              -1.5 * y_per_inch,
              gene.starts[i,'Group.1'], 'black', cex.text=1.2, height=hbar.height)    
  }
  
  par(xpd=F)
  
  
  # draw exons
  #rinfo = cov.info[cov.info$gene == gene,]
  for(e in unique(region.info$exon.start)) {
    abline(v=match(e,region.info$exon.start), col="#aaaaaa")
  }
  
    
  # Add any variants provided to the plot
  if((typeof(variants)=="list") && (nrow(variants)>0)) {
      #print(sprintf("Adding %d variants to plot", nrow(variants)))
      variants.x = sapply(variants$pos, pos_to_x)
      
      # Scale to the maximum depth
      max_depth = max(variants$dp)
      max_bar_height = y.max / 5
      bar_width = 0.1 * x_per_inch 
      bar_top = y.max * 1.035
      alt_height = max_bar_height / max_depth * variants$alt.dp
      rect(
           xleft=variants.x, xright=variants.x+bar_width,
           ytop=rep(bar_top,nrow(variants)), ybottom=rep(bar_top - alt_height, nrow(variants)),
           col="green"
          )
      rect(
           xleft=variants.x, xright=variants.x+bar_width,
           ytop=bar_top-alt_height,
           ybottom=bar_top - alt_height - max_bar_height / max_depth * variants$ref.dp,
           col="red"
          )
  }

  if(!missing(highlight.amplicons)) {
    sample.amplicons = subset(highlight.amplicons, id==s)
    sample.amplicons.overlaps = queryHits(findOverlaps(sample.amplicons, cnvs))
    for(i in sample.amplicons.overlaps) {
      amplicon = sample.amplicons[i]
      plot.hbar(x.pos=c(pos_to_x(start(amplicon)), pos_to_x(end(amplicon))), y.pos=jitter(.25), label="")
    }
  }
}

# This creates a very similar plot to the version above,
# but it shows amplicons that are significantely different in
# coverage drawn in blue. It depends on some 
# external data being available in the environment such as
# 
#   dsd.cnv.candidates
#   dsd.normcov
#   dsd.cov.info
#   dsd.samples
#
plot.cnv.region.amplicons = function(s, gene, cnv.candidates, cnv.cov, dsd.samples, ...) {
  range = cnv.candidates[[s]]
  smoothing = 0.03
  
  cov.info = cnv.cov$cov.info
  cov.norm = cnv.cov$cov.norm
  
  print(fmt("Gene is <%=gene%>"))
  
  # Pull out whole gene
  region = cov.norm[cov.info$gene == gene,]
  rinfo = cov.info[cov.info$gene == gene,]
  
  # Draw the general coverage over the gene
  plot(lowess(region[,s],f=smoothing), t="l", ylim=c(0,2.0), 
       xlab="Position in Region", ylab="Relative Normalised Coverage",
       lwd=1.5,
       main=fmt("Sample <%=s%> Gene <%=gene%>"), ...)
  
  lines(lowess(1+apply(region,1,sd), f=smoothing), col=2, lty=5)
  lines(lowess(1-apply(region,1,sd), f=smoothing), col=2, lty=5)
  abline(h=1.0, col=2)
  
  # Now draw each CNV for that gene
  print(s)
  for(i in 1:length(cnv.candidates[[s]][["cnvs"]])) {
    
    cnv = cnv.candidates[[s]][["cnvs"]][i]
    cnv.gene = mcols(cnv)$gene
    print(cnv.gene)
    if(cnv.gene == gene) {
      
      cnv.start = match(min(rinfo[rinfo$pos>start(cnv),]$pos),rinfo$pos)
      cnv.end= match(max(rinfo[rinfo$pos<end(cnv),]$pos),rinfo$pos)
      
      print(fmt("CNV from <%=cnv.start%> to <%=cnv.end%>"))
      
      # Draw yellow rectangle over region of actual call
      rect(cnv.start-5,-1,cnv.end+5, 2.5, col="#eeee9944")
      
      # Draw lines indicating the start and end of the coverage contigs (sort of representing exons)
      for(e in unique(rinfo$exon.start)) {
        abline(v=match(e,rinfo$exon.start), col="#aaaaaa")
      }
      
      # Draw a line for each amplicon
      cnv.amplicons = dsd.amplicons[[s]][cnv.candidates[[s]][["amplicons"]],] 
      
      cnv.covs = dsd.all[cnv.candidates[[s]][["amplicons"]],dsd.samples] / dsd.all[cnv.candidates[[s]][["amplicons"]],'mean'] 
      
      cnv.amplicon.ranges = GRanges(seqnames=seqnames(cnv), IRanges(start=cnv.amplicons$start, end=cnv.amplicons$end),
                                    relcov=cnv.covs[,s], 
                                    sd=apply(cnv.covs[,dsd.samples[dsd.samples != s]], 1, sd) )
      cnv.amplicon.ranges = cnv.amplicon.ranges[cnv.amplicon.ranges %over% cnv]
      
      print(fmt("Amplicons for cnv <%=i%>: "))
      print(cnv.amplicon.ranges)
      for(j in 1:length(cnv.amplicon.ranges)) {
        print(fmt("Plotting amplicon <%=j%>"))
        cnv.range = cnv.amplicon.ranges[j]
        x.pos = c(match(start(cnv.range),rinfo$pos), match(end(cnv.range)-1,rinfo$pos))
        y.pos = rep(cnv.range$relcov,2)
        print(sprintf("X position = %f, Y position = %f", x.pos, y.pos))
        plot.hbar(x.pos, cnv.range$relcov, label="",col="blue",lwd=2)
        #lines(x.pos, y.pos, col="blue", lwd=2.0)
        
        # Show standard devation for each bar?
        #lines(rep(mean(x.pos),2), c(y.pos[[1]]-cnv.range$sd, y.pos[[1]]+cnv.range$sd), col="blue", lty=2, lwd=2)
      }
    }
  }
}

# Plot each region to a PDF
plot.all.regions = function(regions, cov.norm, cov.info) {
  pdf(fmt("<%=base.dir%>/work/cnv/merged_plots/merged_cnvs.pdf"))
  par(mar=c(10,5,10,5))
  for(i in 1:length(regions)) {
    region = regions[i]
    cnvs = list(xhmm=find.cnv.overlaps(region,dsd.cnv.xhmm),ed=find.cnv.overlaps(region,dsd.cnv.ed), ex=find.cnv.overlaps(region,dsd.cnv.ex))
    plot.cnv.region(region$sample, cnvs, cov.norm, cov.info)
  }
  dev.off()
}

for.each = function(x, fn) {
  for(i in 1:length(x)) {
    fn(x[[i]])
  }
}

as.granges = function(df) {
  # Convert a data.frame to a GRanges object based on chr, start and end columns
  # 
  # Args:
  #   df  A data frame containing chr, start and end columns
  #
  # Result:
  #   a GRanges object representing the genomic ranges in the data frame
  #
  GRanges(df$chr, ranges=IRanges(start=df$start, end=df$end))
}

find.amplicon.overlaps = function(query,amplicons.ranges) {
  query.range = as.granges(query)
  return(subjectHits(findOverlaps(query.range, amplicons.ranges)))
}

find.cnv.variants = function(vcf, cnv) {
  # Find overlapping variants
  #cnv.variants = dsd.variants[(dsd.variants %over% cnv) & 
  #                            (geno(dsd.variants)$GT[,cnv$sample] != "./.") &
  #                            (geno(dsd.variants)$GQ[,cnv$sample]>10)
  #                           ]
  
  # Find overlapping variants
  cnv.variants = vcf[subjectHits(findOverlaps(cnv, vcf))]
  #msg("Found <%=length(cnv.variants)%> over region")
  
  # Find those that have a non-ref genotype in the sample of interest
  x = cnv.variants[geno(cnv.variants)$GT[,cnv$sample]!="./." & geno(cnv.variants)$GQ[,cnv$sample]>10  ]
  
  #cat("Found ",nrow(cnv.variants)," over cnv in sample ", cnv$sample, "\n")
  
  ref.depth = unlist(lapply(geno(x)$AD[,cnv$sample], function(x) if(is.na(x[[1]])) 1 else x[[1]]))
  alt.depth = unlist(lapply(geno(x)$AD[,cnv$sample], function(x) if(length(x) > 1) x[[2]] else 0))
  
  # Find those that have a non-ref genotype in the sample of interest
  data.frame(pos=start(x), 
             gt=geno(x)$GT[,cnv$sample], 
             variant.qual=geno(x)$GQ[,cnv$sample], 
             dp=geno(x)$DP[,cnv$sample],
             ref.dp=ref.depth,
             alt.dp=alt.depth)
}

annotate_variant_counts = function(vcfs, cnvs) {
  if(length(cnvs)==0)
    return(cnvs)
  
  cnvs$het=0
  cnvs$hom=0
  cnvs$bal=0
  for(i in 1:length(cnvs)) {
    print(sprintf("Processing %d / %d",i, length(cnvs)))
    variant.counts = count.variants(vcfs[[cnvs[i]$sample]],cnvs[i])
    cnvs[i]$het = variant.counts[['het']]
    cnvs[i]$hom= variant.counts[['hom']]
    cnvs[i]$bal = variant.counts[['het.ratio']]
  }
  return(cnvs)
}

count.variants = function(vcf,cnv) {
  variants = find.cnv.variants(vcf,cnv)
  ratios = variants$alt.dp/variants$ref.dp
  list(het=sum(variants$gt=="0/1"),hom=sum(variants$gt=="1/1"), het.ratio=mean(ratios[variants$gt=="0/1"]))
}

sim.callers.labels = list(ex="Excavator",ed="ExomeDepth",xhmm="XHMM",ec="ExomeCopy")
sim.callers.colors= list(ex="orange",ed="blue",ang="green",xhmm="red",ec="purple")
sim.callers.symbols= list(ex=21, ed=22,ang=23,xhmm=24,ec=25)
color.complement = list(orange="black", blue="white", green="black", red="black",purple="white")
