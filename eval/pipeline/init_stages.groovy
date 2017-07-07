// Initialization stages - these just set specific output directories for each 
// analysis tool
init_excavator = { branch.excavator_batch=batch_name }
init_xhmm = { branch.xhmm_batch_name=batch_name }
init_exome_depth = {  }
init_cn_mops = { }
init_conifer = { }

create_analysable_target = {
    doc """
          Different CNV callers have various different constraints on what they
          can handle. This stage filters the target regions by a variety of common
          rules to ensure that all the CNV callers included can analyse the data. Also 
          supported is an arbitrary BED file for excluding any regions specified by the user.
        """

    requires target_bed : "The target region BED file",
             sample_names : "The names of samples to analyse"
             
             
    var min_target_size : 30,
        exclude_regions : ''

    def numSamples = sample_names.size()

    from(target_bed) filter('analysable') {

        groovy """
           INCLUDE_CHROMOSOMES = "${chromosomes.join(/,/)}".split(",")

           targetRegion = new BED("$input.bed", withExtra:true).load()

           excludeRegions = "$exclude_regions" ? new BED("$exclude_regions").load() : null

           chrs = targetRegion*.chr.unique()

           chrCounts = targetRegion.countBy { it.chr }

           filteredTargets = targetRegion.grep { 
               (it.chr in INCLUDE_CHROMOSOMES) && 
               (chrCounts[it.chr] > $numSamples) && 
               (it.to - it.from > $min_target_size) && 
               ((excludeRegions == null) || !it.overlaps(excludeRegions))
           } as Regions

           
           filteredTargets.save("$output.bed", sorted:true)
        ""","small"
   }
   branch.analysable_target = output.bed
}


