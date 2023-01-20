
delfin = {
    
    var DELFIN_CONTROL_DIR : false,
        ximmer : null

    requires sample_names : "The names of the samples to analyse"

    var delfin_max_pc_components : 30,
        delfin_cnv_prior : 0.004,
        delfin_deletion_emit_threshold : 3.3,
        delfin_duplication_emit_threshold: 5.5
    
    output.dir = "$batch_name/delfin"
    
    def cov_files = []
    
    if(DELFIN_CONTROL_DIR)
        cov_files = file(DELFIN_CONTROL_DIR)
            .listFiles()
            .findAll { it.name.endsWith('interval_summary') } 

    def test_samples = sample_names
    if(ximmer?.cfg?.containsKey('controls')) {
        test_samples = test_samples.findAll { !(it in ximmer.cfg.controls) }
    }
            
    produce(batch_name + '.delfin.cnvs.tsv') {
        exec """
            $JAVA -Xmx${memory}g -cp '$GROOVY_HOME/lib/*:$GNGS_JAR' gngs.tools.Delfin
                -t $input.bed ${test_samples.collect { "-s $it"}.join(" ")}
                -o $output.cnvs.tsv
                -maxpc $delfin_max_pc_components
                -prior $delfin_cnv_prior
                -del_lrt $delfin_deletion_emit_threshold
                -dup_lrt $delfin_duplication_emit_threshold
                -lr $output.cnvs.tsv.prefix ${cov_files.collect { "-cov $it"}.join(' ')}
                -cov $input.sample_interval_summary
        """, "delfin"
    }
    
    
    branch.caller_result = output.tsv
}
