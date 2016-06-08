// vim: ts=4 sw=4 expandtab

conifer_rpkm = {

    output.dir="$output.dir/rpkms"


    println "Input.bam = $input.bam"


    def rpkmOutput = new SAM(input.bam).samples[0]+".rpkm"

    println "Output rpkm = $rpkmOutput"
    
    produce(rpkmOutput) {
        exec """
        LD_LIBRARY_PATH=$HDF5_DIR/lib $PYTHON $CONIFER/conifer.py rpkm 
          --probes $input.bed
          --output $output.rpkm
          --input $input.bam
        
        ""","conifer"
    }
}

conifer_analyze = {

    requires batch_name : "The name of the batch that Conifer is analysing"

    produce(batch_name+".conifer.hdf5",batch_name+".scree.png") {
        exec """

        LD_LIBRARY_PATH=$HDF5_DIR/lib $PYTHON $CONIFER/conifer.py analyze 
          --probes $input.bed
          --rpkm_dir ${file(input.rpkm).parentFile.absolutePath}
          --output $output.hdf5
          --svd 1
          --write_svals singular_values.txt
          --write_sd sd_values.txt
          --plot_scree $output.png
        """
    }
}

conifer_call = {
    requires batch_name : "The name of the batch that Conifer is analysing"
    produce(batch_name+ ".conifer.cnvs.tsv") {
        exec """
        LD_LIBRARY_PATH=$HDF5_DIR/lib $PYTHON $CONIFER call 
        --input $input.hdf5 
        --output $output.tsv
        """
   }
}

run_conifer = segment {
    create_conifer_target + '%.bam' * [ conifer_rpkm ] + conifer_analyze + conifer_call
}
