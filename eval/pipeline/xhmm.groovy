// vim: ts=4:expandtab:sw=4:cindent
//////////////////////////////////////////////////////////////////
// 
// XHMM Support Routines
//
//////////////////////////////////////////////////////////////////
//
// An example of a pipeline that can be constructed from these:
//
//
//////////////////////////////////////////////////////////////////

// NOTE: Depends on base HG19

xhmm_init = {
    

    // Settings that can be overridden on command line
    var exome_wide_cnv_rate :  '1e-03',
        mean_number_of_targets_in_cnv : 6 ,
        mean_distance_between_targets_within_cnv  : 70,
        mean_of_deletion_z_score_distribution : -3,
        standard_deviation_of_deletion_z_score_distribution: 1,
        mean_of_diploid_z_score_distribution: 0 ,
        standard_deviation_of_diploid_z_score_distribution: 1,
        mean_of_duplication_z_score_distribution: 3,
        standard_deviation_of_duplication_z_score_distribution: 1,
        min_sample_mean:15


    var xhmm_batch_name : batch_name

    branch.xhmm_params = xhmm_batch_name + '.params.txt'
    branch.min_sample_mean = min_sample_mean

    println "Using XHMM ($caller_label) exome wide CNV rate = $exome_wide_cnv_rate Output dir="+output.dir

    produce(xhmm_params) {
          exec """
              printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s" 
                $exome_wide_cnv_rate
                $mean_number_of_targets_in_cnv
                $mean_distance_between_targets_within_cnv
                $mean_of_deletion_z_score_distribution
                $standard_deviation_of_deletion_z_score_distribution
                $mean_of_diploid_z_score_distribution
                $standard_deviation_of_diploid_z_score_distribution
                $mean_of_duplication_z_score_distribution
                $standard_deviation_of_duplication_z_score_distribution
                > $output.txt
          ""","local"
    }
}

gatk_depth_of_coverage = {
    
    // requires target_bed : "BED file containing regions to calculate coverage depth for"

    output.dir = "common/xhmm"
    
    transform("sample_interval_summary") {
        exec """
            $JAVA -Xmx2g -Djava.io.tmpdir=$TMPDIR -jar $GATK/GenomeAnalysisTK.jar 
                 -T DepthOfCoverage 
                 ${inputs.bam.withFlag("-I")}
                 -L $input.bed
                 -R $HGFA
                 -dt BY_SAMPLE 
                 -dcov 5000 
                 -l INFO 
                 --omitDepthOutputAtEachBase 
                 --omitLocusTable 
                 --minBaseQuality 0 
                 --minMappingQuality 20 
                 --start 1 
                 --stop 5000 
                 --nBins 200 
                 --includeRefNSites 
                 --countType COUNT_FRAGMENTS 
                 -o $output.prefix
        """, "gatk_doc_smp"
    }
}

find_extreme_gc_content = {

    // requires target_bed : "BED file containing regions to calculate coverage depth for"

    produce(file(input.bed).name + ".gc.txt", file(input.bed).name+".extremegc.txt") {

        exec """
            $JAVA -Xmx3g -jar $GATK/GenomeAnalysisTK.jar 
                -T GCContentByInterval 
                -L $input.bed
                -R $HGFA
                -o $output1.txt

            cat $output1.txt | awk '{if (\$2 < 0.1 || \$2 > 0.9) print \$1}' > $output2.txt
        ""","medium"
    }
}

xhmm_mean_center = {

    var xhmm_max_mean_rd : 1000

    from("sample_interval_summary", "extremegc.txt") {
        transform("centered","excluded.targets","excluded.samples") {
            exec """
                $XHMM --matrix -r $input.sample_interval_summary
                    --centerData 
                    --centerType target 
                    -o $output.centered
                    --outputExcludedTargets $output.targets
                    --outputExcludedSamples $output.samples
                    --excludeTargets $input2 
                    --minTargetSize 10
                    --maxTargetSize 10000 
                    --minMeanTargetRD 10 --maxMeanTargetRD $xhmm_max_mean_rd 
                    --minMeanSampleRD $min_sample_mean --maxMeanSampleRD $xhmm_max_mean_rd 
                    --maxSdSampleRD 180
            ""","medium"
        }
    }
}

xhmm_merge_coverage = {
    from("*.sample_interval_summary") {
        filter("merged") {
            exec """
                $XHMM --mergeGATKdepths -o $output.sample_interval_summary ${inputs.sample_interval_summary.withFlag("--GATKdepths ")}
            ""","smallish"
        }
    }
}

xhmm_pca = {

    doc "Runs PCA on mean-centered data"

    transform("PC.txt", "PC_SD.txt", "PC_LOADINGS.txt") {
        exec """
            $XHMM --PCA -r $input.centered --PCAfiles $output.txt.prefix.prefix
        """
    }
}

xhmm_normalize = {

    doc "Normalizes mean-centered data using PCA information"
    
    var xhmm_pve_mean_factor : "0.7"

    filter("norm") {
        from("PC.txt") {
            exec """
                $XHMM --normalize 
                   -r $input.centered 
                   --PCAfiles $input.txt.prefix.prefix
                   --normalizeOutput $output.txt
                   --PCnormalizeMethod PVE_mean 
                   --PVE_mean_factor $xhmm_pve_mean_factor
            """
        }
    }
}

xhmm_filter_normalized = {
    doc "Filters and z-score centers (by sample) the PCA-normalized data"
    
    var max_sd_target_rd : 30

    transform("zscored","excluded.targets","excluded.samples") {
        exec """
            $XHMM --matrix -r $input.txt --centerData --centerType sample --zScoreData 
                -o $output.zscored                
                --outputExcludedTargets $output.targets
                --outputExcludedSamples $output.samples
                --maxSdTargetRD $max_sd_target_rd
        """
    }
}

xhmm_filter_orig = {
    doc "Filters original read-depth data to be the same as filtered, normalized data"
    from( "merged.sample_interval_summary", 
          "merged.excluded.targets", 
          "norm.excluded.targets", 
          "merged.excluded.samples", 
          "norm.excluded.samples") {
        filter("filt") {
            exec """
                $XHMM --matrix -r $input.sample_interval_summary
                --excludeTargets $input2
                --excludeTargets $input3
                --excludeSamples $input4
                --excludeSamples $input5
                -o $output.sample_interval_summary
            """
        }

        /*
            --excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt 
            --excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt 
            --excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt 
            --excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt 
        */

    }
}


xhmm_discover = {

    doc "Discovers CNVs in normalized data"

    // requires xhmm_params : "XHMM parameters file, usually 'params.txt'"

    from(xhmm_params) {
        exec """
            $XHMM --discover 
                      -p $input.txt
                      -r $input.zscored 
                      -R $input.sample_interval_summary
                      -c $output.xcnv 
                      -a $output.aux_xcnv 
        """
    }
    
    branch.caller_result = output.xcnv
}

xhmm_pipeline = segment {
    find_extreme_gc_content + 
             xhmm_init + 
             xhmm_merge_coverage + 
             xhmm_mean_center +
             xhmm_pca +
             xhmm_normalize + 
             xhmm_filter_normalized + 
             xhmm_filter_orig + 
             xhmm_discover
}


