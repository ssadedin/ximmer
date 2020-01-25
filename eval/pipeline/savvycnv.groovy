// 
// Implementation of SavvyCNV pipeline
//
init_savvycnv = {
    branch.dir = "analysis/savvy"
}

savvy_bin_coverage = {
    exec """
        export CLASSPATH="$GATK_JAR:$SAVVYCNV_HOME"

        $JAVA -Xmx1g CoverageBinner $input.bam > $output.coverageBinner
    """
}

savvy_call_cnvs = {
    
    var savvy_chunk_size : 200000
    
    exec """
        export CLASSPATH="$GATK_JAR:$SAVVYCNV_HOME"

        $JAVA -Xmx30g SavvyCNV -d $savvy_chunk_size $inputs.coverageBinner > $output.csv 
    """
}

savvy_cnv = segment {
    init_savvycnv + '%.bam' * [ savvy_bin_coverage ] + savvy_call_cnvs
}







