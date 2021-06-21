
delfin = {
    
    requires DELFIN_CONTROL_DIR : "Directory containing interval files of samples to use as controls",
             sample_names : "The names of the samples to analyse"

    var delfin_max_pc_components : 30
    
    output.dir = "analysis/delfin"
    
    def cov_files = 
        file(DELFIN_CONTROL_DIR)
            .listFiles()
            .findAll { it.name.endsWith('interval_summary') } 

    def test_samples = sample_names
    if(ximmer.cfg.containsKey('controls')) {
        test_samples = test_samples.findAll { !(it in ximmer.cfg.controls) }
    }
            
    produce(batch_name + '.delfin.cnvs.tsv') {
        exec """
            $JAVA -Xmx${memory}g -cp $GROOVY_ALL_JAR:$GNGS_JAR gngs.tools.Delfin
                -t $input.bed ${test_samples.collect { "-s $it"}.join(" ")}
                -o $output.cnvs.tsv
                -maxpc $delfin_max_pc_components
                -lr $output.cnvs.tsv.prefix ${cov_files.collect { "-cov $it"}.join(' ')}
                -cov $input.sample_interval_summary
        """, "delfin"
    }
    
    
    branch.caller_result = output.tsv
}
