// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// ExomeDepth Support Routines
//
//////////////////////////////////////////////////////////////////

// By default, run ExomeDepth on each chromosome separately
// This prevents any issues with single-ploidy X chromosomes in males
// but means that chromosomes with small numbers of targets and large events
// could be missed


run_exome_depth = {

    requires target_bed : "BED file containing regions to analyse",
            filtered_sample_names : "List of sample names to process (comma separated, or List object)"

    println "Exome Depth sample names: $sample_names"
    println "Exome Depth branch sample names: $branch.sample_names"
    println "Exome Depth filtered sample names: $filtered_sample_names"

    var transition_probability : "0.0001",
        expected_cnv_length: 50000,
        filter_target_bed : true,
        exome_depth_split_chrs : true,
        filter_to_sex : false
        
    def target_region_to_use = analysable_target
    if(!filter_target_bed) {
        target_region_to_use = target_bed
    }

    def chr = branch.name

    if(filter_to_sex == "FEMALE" && chr == "Y" || chr == "chrY") {
        println "Ignoring Y chromosome due to sex filtering to female samples"
        return
    }
    else {
        println "Analysing $chr with ExomeDepth"
    }

    
    def sample_list = filtered_sample_names
    if(filtered_sample_names instanceof String) {
        sample_list = filtered_sample_names.split(",")
    }

    println "Using Exome Depth transition probability = $transition_probability"

    List outputFiles
    if(exome_depth_split_chrs) {
       outputFiles = [batch_name + '.' + chr + '.exome_depth.tsv', batch_name + '.' + chr + '.exome_depth.warnings.tsv']
    }
    else {
       outputFiles = [batch_name + '.exome_depth.tsv', batch_name + '.exome_depth.warnings.tsv']
    }
    
    produce(outputFiles) {
        R({"""

            source("$TOOLS/r-utils/cnv_utils.R")

            library(ExomeDepth)

            # Reference sequence
            hg19.fasta = "$HGFA"

            # Read the target / covered region
            print(sprintf("Reading target regions for $chr from $target_region_to_use"))

            target.regions = read.bed.ranges(pipe("${exome_depth_split_chrs?"grep '^$chr[^0-9]' $target_region_to_use" : "cat $target_region_to_use"}"))

            # Overlapping targets cause incorrect calls due to ordering applied inside the read counting functions
            # To avoid that, we flatten the target regions here
            targets.flattened = reduce(target.regions)

            # ExomeDepth wants the columns named in a specific way
            dsd.covered = data.frame(
              chromosome=seqnames(targets.flattened),
              start=start(targets.flattened),
              end=end(targets.flattened),
              name=paste(seqnames(targets.flattened),start(targets.flattened),end(targets.flattened),sep="-")
            )

            # Now we need all the bam files. Generate them from sample names
            dsd.samples = c(${sample_list.collect{'"'+it+'"'}.join(",")})

            print(sprintf("Read %d samples",length(dsd.samples)))

            # Here we rely on ASSUMPTIONs:  - Single BAM file per sample
            dsd.bam.files = c(${sample_list.collect { s -> "'$s'='${sample_info[s].files.bam[0]}'"}.join(",") })

            print(sprintf("Found %d bam files",length(dsd.bam.files)))

            # Finally we can call ExomeDepth
            dsd.counts <- getBamCounts(bed.frame = dsd.covered,
                                      bam.files = dsd.bam.files,
                                      include.chr = F,
                                      referenceFasta = hg19.fasta)

            # Old versions of ExomeDepth return IRanges here, newer versions a pure data frame
            # To be more flexible, convert to data frame here. 
            dsd.counts = as.data.frame(dsd.counts)

            # Note: at this point dsd.counts has column names reflecting the file names => convert to actual sample names
            print(sprintf("Successfully counted reads in BAM files"))

            non.sample.columns = (length(colnames(dsd.counts)) - length(dsd.samples))

            colnames(dsd.counts) = c(colnames(dsd.counts)[1:non.sample.columns], dsd.samples)

            # Problem: sample names starting with numbers get mangled. So convert them back, but ignore the first column
            # which is actually the GC percentage
            dsd.samples = colnames(dsd.counts)[-1:-non.sample.columns]

            write(paste("start.p","end.p","type","nexons","start","end","chromosome","id","BF","reads.expected","reads.observed","reads.ratio","sample",sep="\\t"), "$output.exome_depth.tsv")

            for(dsd.test.sample in dsd.samples) {

                print(sprintf("Processing sample %s", dsd.test.sample))

                dsd.reference.samples = dsd.samples[-match(dsd.test.sample, dsd.samples)]

                dsd.counts.df = as.data.frame(dsd.counts[,dsd.reference.samples])

                dsd.test.sample.counts = dsd.counts[,dsd.test.sample]

                assign("last.warning", NULL, envir = baseenv())

                print(sprintf("Selecting reference set for %s ...", dsd.test.sample ))
                dsd.reference = select.reference.set(
                                         test.counts = dsd.test.sample.counts,
                                         reference.counts = as.matrix(dsd.counts.df),
                                         bin.length = dsd.covered$end - dsd.covered$start
                                        )

                all_warnings = data.frame(sample=c(), warning=c())
                if(length(warnings()) > 0) {
                    all_warnings = rbind(all_warnings, data.frame(sample=dsd.test.sample, warning=paste0(warnings(), ',', collapse = '')))
                }
 
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
                                        name = dsd.covered$name,
                                        expected.CNV.length=$expected_cnv_length)


                dsd.results = dsd.cnvs@CNV.calls
                dsd.results$sample = rep(dsd.test.sample, nrow(dsd.results))

                print(sprintf("Writing results ..."))
                if(nrow(dsd.results)>0) {
                    write.table(file="$output.exome_depth.tsv", 
                                x=dsd.results,
                                row.names=F,
                                col.names=F,
                                sep="\\t",
                                append=T)
                } 

                if(nrow(all_warnings)>0) {
                    message(sprintf("Writing %d warnings to warnings file", nrow(all_warnings)))
                }

                write.table(file="$output.warnings.tsv", 
                            x=all_warnings,
                            row.names=F,
                            col.names=F,
                            sep="\t",
                            append=T)    
            }

            print(sprintf("Finished"))

        """},'exome_depth')
    }
    
    if(!exome_depth_split_chrs) {
        branch.caller_result = output.exome_depth.tsv
    }
}

merge_ed = {
    exec """
        cat $inputs.exome_depth.tsv | grep -v '^"sample"' | awk '{ if(NR==1 || \$1 != "start.p") print \$0 }' > ${output(batch_name + ".exome_depth.cnvs.tsv")}
    """, "local"
    
    branch.caller_result = output.tsv
}

exome_depth_pipeline = segment {
    exome_depth_split_chrs ? (chromosomes * [ run_exome_depth ] + merge_ed) : run_exome_depth
}

