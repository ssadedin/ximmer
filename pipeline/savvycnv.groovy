// 
// Implementation of SavvyCNV pipeline
//

init_savvycnv = {
    
    def CONFIG_FILE = "pipeline/config.groovy"
    
    requires SAVVYCNV_HOME : 'The path to SavvyCNV. If you are running this pipeline from Ximmer, ensure the location is configured in $CONFIG_FILE'
    
    if(!file(SAVVYCNV_HOME).exists())
        throw new bpipe.PipelineError("SavvyCNV location $SAVVYCNV_HOME was not found.  Please check this in $CONFIG_FILE")

    branch.SAVVYCNV_JAR="$SAVVYCNV_HOME/build/libs/SavvySuite-all.jar"

    if(!file(SAVVYCNV_JAR).exists())
        throw new bpipe.PipelineError("SavvyCNV jar file $SAVVYCNV_JAR was not found.  Please check SavvySuite was built successfully")

    branch.dir = "analysis/savvy"
}

savvy_bin_coverage = {

    def sample = new gngs.SAM(input.bam.toString()).samples[0]

    produce(sample + '.coverageBinner') {
        exec """
            export CLASSPATH="$SAVVYCNV_JAR"

            $JAVA -Xmx1g CoverageBinner $input.bam > $output.coverageBinner
        ""","savvycnv_coverage"
    }
}

savvy_call_cnvs = {
    
    var savvy_chunk_size : 5000,
        savvy_transition_probability : 0.001
    
    produce('savvy.cnvs.tsv') {
        exec """
            export CLASSPATH="$SAVVYCNV_JAR"

            printf "chromosome\\tstart\\tend\\ttype\\tblock_support\\tblock_span\\t\\tqual\\tqual_rel\\tsample\\n" > $output.tsv

            $JAVA -Xmx${memory}g SavvyCNV -data -d $savvy_chunk_size -trans $savvy_transition_probability $inputs.coverageBinner >> $output.tsv 
        ""","savvycnv"
    }

    branch.caller_result = output.tsv
}

savvy_cnv = segment {
    init_savvycnv + '%.bam' * [ savvy_bin_coverage ] + savvy_call_cnvs
}
