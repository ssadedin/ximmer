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


exome_depth_count_fragments = {

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
    
    def outputFile = exome_depth_split_chrs ? 
        batch_name + '.' + chr + '.counts.tsv.gz' 
       :
        batch_name + '.counts.tsv.gz' 
       

    produce(batch_name + '.' + chr + '.counts.tsv.gz') {
        R({"""

            source("$TOOLS/r-utils/cnv_utils.R")

            library(ExomeDepth)

            # Reference sequence
            ref.fasta = "$HGFA"

            # Read the target / covered region
            print(sprintf("Reading target regions for $chr from $target_region_to_use"))

            target.regions = read.bed.ranges(pipe("${exome_depth_split_chrs?"grep '^$chr[^0-9]' $target_region_to_use" : "cat $target_region_to_use"}"))

            # Overlapping targets cause incorrect calls due to ordering applied inside the read counting functions
            # To avoid that, we flatten the target regions here
            targets.flattened = reduce(target.regions)

            # ExomeDepth wants the columns named in a specific way
            target.covered = data.frame(
              chromosome=seqnames(targets.flattened),
              start=start(targets.flattened),
              end=end(targets.flattened),
              name=paste(seqnames(targets.flattened),start(targets.flattened),end(targets.flattened),sep="-")
            )

            # Now we need all the bam files. Generate them from sample names
            ed.samples = c(${sample_list.collect{'"'+it+'"'}.join(",")})

            print(sprintf("Read %d samples",length(ed.samples)))

            # Here we rely on ASSUMPTIONs:  - Single BAM file per sample
            ed.bam.files = c(${sample_list.collect { s -> "'$s'='${sample_info[s].files.bam[0]}'"}.join(",") })

            print(sprintf("Found %d bam files",length(ed.bam.files)))

            # Finally we can call ExomeDepth
            ed.counts <- getBamCounts(bed.frame = target.covered,
                                      bam.files = ed.bam.files,
                                      include.chr = F,
                                      referenceFasta = ref.fasta)

            # Old versions of ExomeDepth return IRanges here, newer versions a pure data frame
            # To be more flexible, convert to data frame here. 
            ed.counts = as.data.frame(ed.counts)

            # Note: at this point ed.counts has column names reflecting the file names => convert to actual sample names
            print(sprintf("Successfully counted reads in BAM files"))

            non.sample.columns = (length(colnames(ed.counts)) - length(ed.samples))

            colnames(ed.counts) = c(colnames(ed.counts)[1:non.sample.columns], ed.samples)

            count.file = gzfile('$output.gz','w')
            write.table(ed.counts, file=count.file, row.names=FALSE)
            close(count.file)
        """},'exome_depth')
    }
}

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
    
    List chrPart = exome_depth_split_chrs ? [chr] : []
    
    List outputFiles = [
       [ batch_name ] + chrPart + ['exome_depth.tsv'],
       [ batch_name ] + chrPart + ['exome_depth.warnings.tsv'],
       [ batch_name ] + chrPart + ['exome_depth.refstats.tsv']      
    ]*.join('.')
   
    produce(outputFiles) {
        R({"""
            source("$TOOLS/r-utils/cnv_utils.R")

            library(ExomeDepth)

            # Reference sequence
            ref.fasta = "$HGFA"

            # Read the target / covered region
            print(sprintf("Reading target regions for $chr from $target_region_to_use"))

            target.regions = read.bed.ranges(pipe("${exome_depth_split_chrs?"grep '^$chr[^0-9]' $target_region_to_use" : "cat $target_region_to_use"}"))

            # Overlapping targets cause incorrect calls due to ordering applied inside the read counting functions
            # To avoid that, we flatten the target regions here
            targets.flattened = reduce(target.regions)

            # ExomeDepth wants the columns named in a specific way
            target.covered = data.frame(
              chromosome=seqnames(targets.flattened),
              start=start(targets.flattened),
              end=end(targets.flattened),
              name=paste(seqnames(targets.flattened),start(targets.flattened),end(targets.flattened),sep="-")
            )

            # Now we need all the bam files. Generate them from sample names
            ed.samples = c(${sample_list.collect{'"'+it+'"'}.join(",")})
            ed.test.samples = c(${test_samples.collect{'"'+it+'"'}.join(",")})

            print(sprintf("Read %d samples",length(ed.samples)))


            count.file = gzfile('$input.counts.tsv.gz','r')
            ed.counts = read.table(count.file, header=TRUE, stringsAsFactors=FALSE)
            close(count.file)

            print("Read counts from $input.gz for $chr")

            non.sample.columns = ncol(ed.counts) - length(ed.samples)

            colnames(ed.counts) = c(names(ed.counts)[1:non.sample.columns], ed.samples)

            write(paste("start.p","end.p","type","nexons","start","end","chromosome","id","BF","reads.expected","reads.observed","reads.ratio","sample",sep="\\t"), "$output.exome_depth.tsv")

            reference.choices = data.frame(list(sample=c(), choices=c()))

            all.reference.stats = NA

            for(ed.test.sample in ed.test.samples) {

                print(sprintf("Processing sample %s", ed.test.sample))

                reference.set.samples = ed.samples[-match(ed.test.sample, ed.samples)]

                sample.reference.counts = as.data.frame(ed.counts[,reference.set.samples])

                ed.test.sample.counts = ed.counts[,ed.test.sample]

                #assign("last.warning", NULL, envir = baseenv())

                print(sprintf("Selecting reference set for %s ...", ed.test.sample ))
                reference.set = select.reference.set(
                                         test.counts = ed.test.sample.counts,
                                         reference.counts = as.matrix(sample.reference.counts),
                                         bin.length = target.covered\$end - target.covered\$start
                                        )

                sample.reference.stats = as.data.frame(reference.set\$summary.stats)
                sample.reference.stats$sample = ed.test.sample 
             
                if(is.na(all.reference.stats)) {
                    all.reference.stats = sample.reference.stats
                }
                else {
                    all.reference.stats = rbind(all.reference.stats, sample.reference.stats)
                }

                all_warnings = data.frame(sample=c(), warning=c())
                if(length(warnings()) > 0) {
                    all_warnings = rbind(all_warnings, data.frame(sample=ed.test.sample, warning=paste0(warnings(), ',', collapse = '')))
                }
 
                # Get counts just for the reference set
                reference.set.counts = apply(sample.reference.counts[,reference.set\$reference.choice,drop=F],1,sum)

                print(sprintf("Creating ExomeDepth object ..."))
                sample.ed = new("ExomeDepth",
                              test = ed.test.sample.counts,
                              reference = reference.set.counts,
                              formula = "cbind(test, reference) ~ 1")


                print(sprintf("Calling CNVs ..."))
                sample.cnvs = CallCNVs(x = sample.ed,
                                        transition.probability = $transition_probability,
                                        chromosome = target.covered\$chromosome,
                                        start = target.covered\$start,
                                        end = target.covered\$end,
                                        name = target.covered\$name,
                                        expected.CNV.length=$expected_cnv_length)


                sample.results = sample.cnvs@CNV.calls
                sample.results$sample = rep(ed.test.sample, nrow(sample.results))

                print(sprintf("Writing results ..."))
                if(nrow(sample.results)>0) {
                    write.table(file="$output.exome_depth.tsv", x=sample.results, row.names=F, col.names=F, sep="\\t", append=T)
                } 

                if(nrow(all_warnings)>0) {
                    message(sprintf("Writing %d warnings to warnings file", nrow(all_warnings)))
                }

                write.table(file="$output.warnings.tsv", x=all_warnings, row.names=F, col.names=F, sep="\t", append=T)    
            }

            write.table(file="$output.refstats.tsv", x=all.reference.stats, row.names=F, col.names=T, sep="\t")    

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
    if(exome_depth_split_chrs)
        return (chromosomes * [ exome_depth_count_fragments + run_exome_depth ] + merge_ed) 
    else
        return (exome_depth_count_fragments + run_exome_depth)
}

