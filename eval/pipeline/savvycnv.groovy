// 
// Implementation of SavvyCNV pipeline
//
init_savvycnv = {
    branch.dir = "analysis/savvy"
}

savvy_bin_coverage = {

    requires GATK4_JAR : 'Location of jar file for GATK'

    def sample = new gngs.SAM(input.bam.toString()).samples[0]

    produce(sample + '.coverageBinner') {
        exec """
            export CLASSPATH="$GATK4_JAR:$SAVVYCNV_HOME"

            $JAVA -Xmx1g CoverageBinner $input.bam > $output.coverageBinner
        """
    }
}

savvy_call_cnvs = {
    
    var savvy_chunk_size : 30000,
        savvy_transition_probability : 0.01
    
    produce('savvy.cnvs.tsv') {
        exec """
            export CLASSPATH="$GATK4_JAR:$SAVVYCNV_HOME"

            printf "chromosome\\tstart\\tend\\ttype\\tblock_support\\tblock_span\\t\\tqual\\tqual_rel\\tsample\\n" > $output.tsv

            $JAVA -Xmx30g SavvyCNV -data -d $savvy_chunk_size -trans $savvy_transition_probability $inputs.coverageBinner >> $output.tsv 
        """
    }

    branch.caller_result = output.tsv
}

savvy_cnv = segment {
    init_savvycnv + '%.bam' * [ savvy_bin_coverage ] + savvy_call_cnvs
}







