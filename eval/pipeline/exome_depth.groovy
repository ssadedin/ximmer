// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// ExomeDepth Support Routines
//
//////////////////////////////////////////////////////////////////
run_exome_depth = {

    requires target_bed : "BED file containing regions to analyse",
            sample_names : "List of sample names to process (comma separated, or List object)",
            chr : "Chromosome to process for"

    var transition_probability : "0.0001"

    def chr = branch.name
    
    def sample_list = sample_names
    if(sample_names instanceof String) {
        sample_list = sample_names.split(",")
    }

    println "Using Exome Depth transition probability = $transition_probability"

    R({"""

        source("$TOOLS/r-utils/cnv_utils.R")

        library(ExomeDepth)

        # Reference sequence
        hg19.fasta = "$HGFA"

        # Read the target / covered region
        print(sprintf("Reading target regions for $chr from $target_bed"))
        dsd.covered = read.bed(pipe("grep '$chr[^0-9]' $target_bed"))

        # ExomeDepth wants the columns named differently
        dsd.covered = data.frame(
            chromosome=dsd.covered\$chr, 
            start=dsd.covered$start, 
            end=dsd.covered$end,  
            name=paste(dsd.covered\$chr,dsd.covered$start,dsd.covered\$end,sep="-")
        )

        # Now we need all the bam files. Generate them from sample names
        dsd.samples = c(${sample_list.collect{'"'+it+'"'}.join(",")})

        print(sprintf("Read %d samples",length(dsd.samples)))

        # Here we rely on ASSUMPTIONs:  - Single BAM file per sample
        dsd.bam.files = c(${sample_info.collect { key, s -> "'$s.sample'='${s.files.bam[0]}'"}.join(",") })

        print(sprintf("Found %d bam files",length(dsd.bam.files)))

        # Finally we can call ExomeDepth
        dsd.counts <- getBamCounts(bed.frame = dsd.covered,
                                  bam.files = dsd.bam.files,
                                  include.chr = F,
                                  referenceFasta = hg19.fasta)

        # Note: at this point dsd.counts has column names reflecting the file names => convert to actual sample names
        print(sprintf("Successfully counted reads in BAM files"))

        colnames(dsd.counts) = c("GC", dsd.samples)

        # Problem: sample names starting with numbers get mangled. So convert them back, but ignore the first column
        # which is actually the GC percentage
        dsd.samples = colnames(dsd.counts)[-1]

	write(paste("start.p","end.p","type","nexons","start","end","chromosome","id","BF","reads.expected","reads.observed","reads.ratio","sample",sep="\\t"), "$output.tsv")

        for(dsd.test.sample in dsd.samples) {

            print(sprintf("Processing sample %s", dsd.test.sample))

            dsd.reference.samples = dsd.samples[-match(dsd.test.sample, dsd.samples)]

            dsd.counts.df = as.data.frame(dsd.counts[,dsd.reference.samples])[,-1:-6]

            dsd.test.sample.counts = dsd.counts[,dsd.test.sample][[1]]

            print(sprintf("Selecting reference set for %s ...", dsd.test.sample ))
            dsd.reference = select.reference.set(
                                     test.counts = dsd.counts[,dsd.test.sample][[1]],
                                     reference.counts = as.matrix(as.data.frame(dsd.counts[,dsd.reference.samples])[,-1:-6]),
                                     bin.length = dsd.covered$end - dsd.covered$start
                                    )

            # Get counts just for the reference set
            dsd.reference.counts = apply(dsd.counts.df[,dsd.reference\$reference.choice,drop=F],1,sum)

            print(sprintf("Creating ExomeDepth object ..."))
            dsd.ed = new("ExomeDepth",
                          test = dsd.test.sample.counts,
                          reference = dsd.reference.counts,
                          formula = "cbind(test, reference) ~ 1")


            print(sprintf("Calling CNVs ..."))
            dsd.cnvs = CallCNVs(x = dsd.ed,
                                    transition.probability = $transition_probability,
                                    chromosome = dsd.covered$chromosome,
                                    start = dsd.covered$start,
                                    end = dsd.covered$end,
                                    name = dsd.covered$name)


            dsd.results = dsd.cnvs@CNV.calls
            dsd.results$sample = rep(dsd.test.sample, nrow(dsd.results))

            print(sprintf("Writing results ..."))
            if(nrow(dsd.results)>0) {
                write.table(file="$output.tsv", 
                            x=dsd.results,
                            row.names=F,
                            col.names=F,
                            append=T)
            } 
        }

        print(sprintf("Finished"))

    """},'exome_depth')
}

merge_ed = {
    exec """
        cat $inputs.exome_depth.tsv | grep -v '^"sample"' | awk '{ if(NR==1 || \$1 != "\\"start.p\\"") print \$0 }' > ${output(batch_name + ".exome_depth.cnvs.tsv")}
    """
}

exome_depth_pipeline = segment {
    chromosomes * [ run_exome_depth ] + merge_ed
}

