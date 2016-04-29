// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// Excavator Support Routines
//
//////////////////////////////////////////////////////////////////
//
// An example of a pipeline that can be constructed from these:
//
//   excavator_target + run_samples * [ create_excavator_config + excavator ] + combine_excavator_results
//
//////////////////////////////////////////////////////////////////

excavator_target = {
    doc "Create metadata about target region based on mappability, GC content and other info"

    requires target_bed : """BED file describing target region on which analysis is to be run. Usually 
                             this should be the whole region covered by probes / amplicons, not just 
                             the actual genes that are of interest, and should be unique"""

   def (mappability, reference) = file("$EXCAVATOR/hg19.sourcetarget.txt").text.split('\\s')
   check {
       exec """
           [ -e $mappability ] && [ -e $reference ]
       """
   } otherwise {
       fail "The file $EXCAVATOR/hg19.sourcetarget.txt contains references to files that don't exist"
   }

    // Compute the target region
    produce("excavator.target.name.txt") {
        exec """
            sha1sum $target_bed | awk '{ print \$1 }' | 
                xargs echo `basename $target_bed` | 
                sed 's/.bed /_/' > $output.txt
        """
    }

    println "Loading excavator target name from $output.txt"

    Thread.sleep(1000)
    File targetNameFile = file(output.txt)
    targetNameFile.parentFile.listFiles()

    if(file(output.txt).exists()) {
        branch.excavator_target_name=file(output.txt).text.trim()
        println "The excavator target name is $excavator_target_name"
    }
    else {
        fail "WARNING: expected target name file $output.txt does not exist!"
    }

    produce("$EXCAVATOR/data/targets/hg19/$excavator_target_name/MAP/Map.chrX.RData") {
        exec """
                cd $EXCAVATOR 

                export PATH=`dirname $SAMTOOLS`:"$PATH"

                perl TargetPerla.pl hg19.sourcetarget.txt $target_bed $excavator_target_name
            """
    }
}

create_excavator_config = {

    doc """
            Build an Excavator configuration file that specifies the current sample
            as a test sample, but the other samples as controls.
        """

    branch.sample = branch.name

    msg "Creating excavator config for $sample now"

    if(!(sample in all_samples))
        fail "Sample '$sample' could not be located in the known samples parsed from the input BAM files"

    produce("dsd.excavator.input.${sample}.txt") {
        if(file(output.toString()).exists())
            return
        new File(output).withWriter { w ->
            msg "Creating $output"

            // First sample as the test candidate
            w.println "${excavator_target_name}\thg19\t${all_samples[sample].files.bam[0]}\t$sample\tT1"
            def count = 1

            // All the other samples as controls
            w.print all_samples.grep { it.key != sample }.collect { 
                "$excavator_target_name\thg19\t${all_samples[it.key].files.bam[0]}\t$it.key\tC${count++}\n" 
            }.join('')
        }
    }
}

run_excavator = {

   var excavator_theta : "1e-4",
       excavator_dnorm : "10e6",
       excavator_min_exons : "2"

    def outputName = output.dir + "/${excavator_batch}_${sample}_outputs/Results/${sample}/FastCallResults_${sample}.txt"
    def exParamFile = output.dir + "/${excavator_batch}_${sample}_outputs/Results/${sample}/params_${sample}.txt"


    def excavatorOutputDir = branch.dir + "/${excavator_batch}_${sample}_outputs"

    // Excavator actually writes the result we want into a subdirectory of the directory you tell it
    output.dir = excavatorOutputDir + "/Results/${sample}"

    // println "Output will be $outputName"

    produce("FastCallResults_${sample}.txt","params_${sample}.txt") {

        file(output.dir).mkdirs()

        file(output2).text = """
            ## Omega parameter for the HSLM algorithm ##
            0.1
            ## Theta parameter (baseline probability m_i changes its value) for the HSLM algorithm ##
            $excavator_theta
            ## D_norm parameter for the HSLM algorithm ##
            $excavator_dnorm
            ## Cellularity parameter for the FastCall Calling algorithm ##
            1
            ## Threshold d for the truncated gaussian distribution of the FastCall Calling algorithm ##
            0.5
            ## Threshold u for the truncated gaussian distribution of the FastCall Calling algorithm ##
            0.35
            ## Segment with a number of exons smaller than a threshold are filtered out ##
            $excavator_min_exons 
        """.stripIndent().trim()

        println "Wrote EXCAVATOR parameters: $output2.txt"

        println  "Creating $output1 and $output2"

        exec """
            export PATH=`dirname $SAMTOOLS`:"$PATH"

            echo "Output is $output1, dir=$output.dir"

            perl $EXCAVATOR/ReadPerla.pl -v $input.txt -p $exParamFile --mode pooling --mapq 1 ${excavatorOutputDir}
        """, "excavator"
    }
}

combine_excavator_results = {

    produce("${excavator_batch}.excavator.cnvs.tsv") {
        exec """
            for i in ${run_samples.join(' ')};
            do 
                echo "Processing sample $i";

                result=${output.dir}/${excavator_batch}_\${i}_outputs/Results/${i}/FastCallResults_${i}.txt ; 

                if [ -e $result ]; 
                then 
                    echo "Found $result"

                    cat $result | grep -v Chromosome | awk '{ print "'$i'" "\\t" \$0 }' >> $output.tsv; 
                fi 
            done
        """
    }
}

excavator_pipeline = segment {
    excavator_target + run_samples * [ create_excavator_config + run_excavator ] + combine_excavator_results
}

