// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// Pipeline stage to run cn.MOPs on exome data
//
//////////////////////////////////////////////////////////////////
cn_mops_call_cnvs = {

    var batch_name : false

    var prior_impact : 10,
        min_width : 5,
        lower_threshold : -0.8

    def outputFile = batch_name ? batch_name + '.cnmops.cnvs.tsv' : input.bam + '.cnmops.cnvs.tsv'

    produce(outputFile) {
        R({"""

            source("$TOOLS/r-utils/cnv_utils.R")

            library(cn.mops)
            library(Rsamtools)

            bam.files = c('${inputs.bam.join("','")}')

            target.region = unique(read.bed.ranges("$input.bed"))

            # We have to order by chromosome lexically, as this is 
            # what countBam does internally - otherwise getSegmentReadCountsFromBAM
            # returns wrongly ordered counts
            target.region = target.region[order(seqnames(target.region))]

            bam.samples = sapply(bam.files, function(file.name) {
                # Scan the bam header and parse out the sample name from the first read group
                read.group.info = strsplit(scanBamHeader(file.name)[[1]]$text[["@RG"]],":")
                names(read.group.info) = sapply(read.group.info, function(field) field[[1]]) 
                return(read.group.info$SM[[2]])
            })


            print(sprintf("Counting reads from %d bam files for %d regions",length(bam.files), length(target.region)))
            mops.counts = getSegmentReadCountsFromBAM(bam.files, GR=target.region,mode="paired")

            print("Normalizing read counts ...")
            mops.counts.norm <- normalizeChromosomes(mops.counts)

            print("Fitting MOPs ...")
            mops.results = exomecn.mops(mops.counts.norm,norm=FALSE,
                    priorImpact=$prior_impact,
                    minWidth=$min_width,
                    lowerThreshold=$lower_threshold)

            mops.results.cn = calcIntegerCopyNumbers(mops.results)

            mops.cnvs = mops.results.cn@cnvs
            mops.cnvs$sampleName = sapply(mops.cnvs$sampleName, function(sn) {
                bam.samples[[sn]]
            })

            print(sprintf("Writing results ..."))
            write.table(file="$output.tsv", 
                            x=as.data.frame(mops.cnvs),
                            row.names=F,
                            quote=F)

        """}, "cnmops")
    }
}
