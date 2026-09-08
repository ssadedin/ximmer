package ximmer

import gngs.Region
import gngs.Regions
import gngs.ToolBase
import gngs.XPos
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.*
import htsjdk.variant.vcf.*

import java.nio.file.Files

/**
 * Merges per-sample CNV/SV VCFs (produced by TSVtoVCF) into a single
 * pseudo-joint-genotyped VCF for a trio (proband + up to 2 parents).
 * 
 * <p>The merge is <b>proband-driven</b>: the proband's event definitions are
 * authoritative. For each proband event, parents are tested for a matching
 * event. In this initial version, matching requires <b>exact</b> boundary
 * agreement (same chr, start, end, and SVTYPE).</p>
 * 
 * <p>Processing order:</p>
 * <ol>
 *   <li>For each proband event, find exact-matching parent events and emit
 *       a joint line. Matched parent events are marked as consumed.</li>
 *   <li>For remaining (unconsumed) parent events, test for exact match
 *       between the two parents. If matched, emit a joint line with proband
 *       as hom-ref. If unmatched, emit a single-parent line with all others
 *       as hom-ref.</li>
 *   <li>All output variants are sorted by genomic position before writing.</li>
 * </ol>
 */
@Log
class FamilyMerge extends ToolBase {

    static void main(String[] args) {
        cli('FamilyMerge -proband <sample> -parent <sample> [-parent <sample2>] -i <vcf> [-i <vcf> ...] -o <output.vcf>',
            'Merges per-sample CNV VCFs into a joint-genotyped family VCF', args) {
            i 'Input VCF file (one per sample, may be specified multiple times)', args: '*', required: true, type: File
            proband 'Proband sample name', args: 1, required: true, type: String
            parent 'Parent sample name (specify once or twice)', args: '*', required: true, type: String
            o 'Output VCF file', args: 1, required: true, type: File
        }
    }

    @Override
    public void run() {
        String probandSample = opts['proband']
        List<String> parentSamples = opts['parents']

        if(parentSamples.size() > 2) {
            throw new IllegalArgumentException("At most 2 parents are supported for trio merging, but ${parentSamples.size()} were specified")
        }

        List<File> inputFiles = opts['is']
        log.info "Merging ${inputFiles.size()} VCF files for proband=${probandSample}, parents=${parentSamples}"

        Map<String, List<VariantContext>> variantsBySample = loadAllVariants(inputFiles)

        List<String> allSamples = [probandSample] + parentSamples

        // Validate that we found variants for the proband
        if(!variantsBySample.containsKey(probandSample)) {
            log.warning("No variants found for proband ${probandSample}")
        }

        // Collect VCF header info from first input
        VCFHeader sourceHeader = readHeader(inputFiles[0])

        opts['o'].withWriter { Writer w ->
            List<VariantContext> merged = mergeVariants(probandSample, parentSamples, variantsBySample)
            writeVCF(w, sourceHeader, allSamples, merged)
        }

        log.info "Wrote ${opts['o']}"
    }

    /**
     * Load all variants from the given VCF files, grouped by sample name.
     * Each VCF is expected to contain genotypes; the non-hom-ref sample is
     * extracted as the "owning" sample for each variant.
     *
     * @param vcfFiles list of input VCF files
     * @return map of sample name to list of VariantContext objects
     */
    Map<String, List<VariantContext>> loadAllVariants(List<File> vcfFiles) {
        Map<String, List<VariantContext>> result = [:]
        for(File f in vcfFiles) {
            log.info "Loading variants from ${f}"
            VCFFileReader reader = new VCFFileReader(f, false)
            try {
                VCFHeader fileHeader = reader.getFileHeader()
                List<String> vcfSamples = fileHeader.getGenotypeSamples()
                if(vcfSamples.isEmpty()) {
                    log.warning("VCF file ${f} has no sample columns, skipping")
                    continue
                }
                // Use the first sample from the header as the owning sample
                String sample = vcfSamples[0]
                log.info "VCF file ${f} contains sample: ${sample}"
                if(!result.containsKey(sample)) {
                    result[sample] = new ArrayList(1000)
                }
                for(VariantContext vc in reader) {
                    result[sample].add(vc)
                }
            }
            finally {
                reader.close()
            }
        }
        for(Map.Entry<String, List<VariantContext>> entry in result) {
            log.info "Loaded ${entry.value.size()} variants for sample ${entry.key}"
        }
        return result
    }

    /**
     * Read the VCF header from a file.
     */
    VCFHeader readHeader(File vcfFile) {
        VCFFileReader reader = new VCFFileReader(vcfFile, false)
        try {
            return reader.getFileHeader()
        }
        finally {
            reader.close()
        }
    }

    /**
     * Find the sample in a VariantContext that has a non-hom-ref genotype.
     *
     * @return sample name, or null if all samples are hom-ref
     */
    @CompileStatic
    static String findNonRefSample(VariantContext vc) {
        for(Genotype g in vc.getGenotypes()) {
            if(!g.isHomRef()) {
                return g.getSampleName()
            }
        }
        return null
    }

    /**
     * Compute the key used for exact-match grouping: (chr, start, end, svtype).
     */
    @CompileStatic
    static List exactMatchKey(VariantContext vc) {
        return [vc.getContig(), vc.getStart(), vc.getEnd(), vc.getAttribute("SVTYPE")]
    }

    /**
     * Core merge logic for a trio.
     *
     * <p>Phase 1: proband-driven. For each proband event, look for exact-matching
     * parent events. Emit a joint line. Mark matched parent events as consumed.</p>
     *
     * <p>Phase 2: remaining parent events. For unconsumed parent events, test
     * for exact match between the two parents. Emit joint or single lines
     * with proband as hom-ref.</p>
     *
     * @param probandSample the proband sample name
     * @param parentSamples list of parent sample names (1 or 2)
     * @param variantsBySample map of sample to variant list
     * @return list of merged VariantContext objects (unsorted)
     */
    List<VariantContext> mergeVariants(String probandSample, List<String> parentSamples,
                                       Map<String, List<VariantContext>> variantsBySample) {

        List<String> allSamples = [probandSample] + parentSamples

        // Index parent variants by exact-match key for O(1) lookup
        // Each entry maps to a list because multiple calls could theoretically share a key
        Map<String, Map<List, List<VariantContext>>> parentIndexes = [:]
        for(String parent in parentSamples) {
            Map<List, List<VariantContext>> idx = [:]
            List<VariantContext> parentVars = variantsBySample[parent] ?: []
            for(VariantContext vc in parentVars) {
                List key = exactMatchKey(vc)
                if(!idx.containsKey(key)) {
                    idx[key] = []
                }
                idx[key].add(vc)
            }
            parentIndexes[parent] = idx
        }

        // Track consumed parent events (by identity)
        Set<VariantContext> consumedParentEvents = new HashSet<VariantContext>()

        List<VariantContext> output = []

        // Phase 1: proband-driven
        List<VariantContext> probandVars = variantsBySample[probandSample] ?: []
        log.info "Phase 1: Processing ${probandVars.size()} proband events"

        for(VariantContext probandVc in probandVars) {
            List key = exactMatchKey(probandVc)

            // Find matching parent events
            Map<String, VariantContext> parentMatches = [:]
            for(String parent in parentSamples) {
                List<VariantContext> candidates = parentIndexes[parent][key]
                if(candidates != null && !candidates.isEmpty()) {
                    // Take the first unconsumed match
                    VariantContext match = candidates.find { !consumedParentEvents.contains(it) }
                    if(match != null) {
                        parentMatches[parent] = match
                        consumedParentEvents.add(match)
                    }
                }
            }

            // Build joint variant
            VariantContext joint = buildJointVariant(probandVc, probandSample, allSamples, parentMatches)
            output.add(joint)
        }

        log.info "Phase 1 complete: ${output.size()} proband-driven variants emitted, ${consumedParentEvents.size()} parent events consumed"

        // Phase 2: remaining parent events
        // Collect unconsumed parent events
        Map<String, List<VariantContext>> unconsumedByParent = [:]
        for(String parent in parentSamples) {
            List<VariantContext> parentVars = variantsBySample[parent] ?: []
            unconsumedByParent[parent] = parentVars.findAll { !consumedParentEvents.contains(it) }
        }

        if(parentSamples.size() == 2) {
            String parent1 = parentSamples[0]
            String parent2 = parentSamples[1]

            // Index unconsumed parent2 events for matching against parent1
            Map<List, List<VariantContext>> parent2Index = [:]
            for(VariantContext vc in unconsumedByParent[parent2]) {
                List key = exactMatchKey(vc)
                if(!parent2Index.containsKey(key)) {
                    parent2Index[key] = []
                }
                parent2Index[key].add(vc)
            }

            Set<VariantContext> consumedParent2 = new HashSet<VariantContext>()

            // Try to match parent1 unconsumed events with parent2
            for(VariantContext p1Vc in unconsumedByParent[parent1]) {
                List key = exactMatchKey(p1Vc)
                List<VariantContext> p2Candidates = parent2Index[key]
                VariantContext p2Match = null
                if(p2Candidates != null && !p2Candidates.isEmpty()) {
                    p2Match = p2Candidates.find { !consumedParent2.contains(it) }
                    if(p2Match != null) {
                        consumedParent2.add(p2Match)
                    }
                }

                Map<String, VariantContext> matches = [(parent1): p1Vc]
                if(p2Match != null) {
                    matches[parent2] = p2Match
                }

                // Primary event is parent1's, proband is hom-ref
                VariantContext joint = buildJointVariant(p1Vc, parent1, allSamples, matches)
                output.add(joint)
            }

            // Emit remaining unmatched parent2 events
            for(VariantContext p2Vc in unconsumedByParent[parent2]) {
                if(!consumedParent2.contains(p2Vc)) {
                    Map<String, VariantContext> matches = [(parent2): p2Vc]
                    VariantContext joint = buildJointVariant(p2Vc, parent2, allSamples, matches)
                    output.add(joint)
                }
            }
        }
        else if(parentSamples.size() == 1) {
            // Single parent - just emit unconsumed events
            String parent = parentSamples[0]
            for(VariantContext pVc in unconsumedByParent[parent]) {
                Map<String, VariantContext> matches = [(parent): pVc]
                VariantContext joint = buildJointVariant(pVc, parent, allSamples, matches)
                output.add(joint)
            }
        }

        int parentOnlyCount = output.size() - probandVars.size()
        log.info "Phase 2 complete: ${parentOnlyCount} parent-only variants emitted"
        log.info "Total output variants: ${output.size()}"

        return output
    }

    /**
     * Build a joint-genotyped VariantContext from a primary event definition
     * and a map of sample matches.
     *
     * <p>The primary event defines the coordinates, alleles, and INFO fields.
     * For each sample in allSamples:</p>
     * <ul>
     *   <li>If the sample has a match in sampleMatches, build a het genotype
     *       with CR and NC from that sample's own call</li>
     *   <li>Otherwise, build a hom-ref genotype</li>
     * </ul>
     *
     * @param primaryVc the variant context that defines the event coordinates
     * @param primarySample the sample that "owns" the primary event
     * @param allSamples all samples to include in the output
     * @param sampleMatches map of sample name to their matching VariantContext (may not include primarySample)
     * @return a new VariantContext with genotypes for all samples
     */
    @CompileStatic
    VariantContext buildJointVariant(VariantContext primaryVc, String primarySample,
                                     List<String> allSamples,
                                     Map<String, VariantContext> sampleMatches) {

        Allele refAllele = primaryVc.getReference()
        Allele altAllele = primaryVc.getAlternateAllele(0)

        List<Genotype> genotypes = new ArrayList(allSamples.size())

        for(String sample in allSamples) {
            VariantContext matchedVc = sampleMatches[sample]
            if(matchedVc != null) {
                // This sample has the event - extract their FORMAT fields
                Genotype sourceGt = matchedVc.getGenotype(sample)
                if(sourceGt != null && !sourceGt.isHomRef()) {
                    // Carry over the original genotype's FORMAT fields
                    Map<String, Object> formatFields = new LinkedHashMap()
                    Object cr = sourceGt.getExtendedAttribute("CR")
                    Object nc = sourceGt.getExtendedAttribute("NC")
                    if(cr != null) {
                        formatFields.put("CR", cr)
                    }
                    if(nc != null) {
                        formatFields.put("NC", nc)
                    }
                    genotypes.add(GenotypeBuilder.create(sample, [refAllele, altAllele], formatFields))
                }
                else {
                    // Sample exists in the matched VC but is hom-ref there (shouldn't normally happen)
                    genotypes.add(GenotypeBuilder.create(sample, [refAllele, refAllele]))
                }
            }
            else if(sample == primarySample) {
                // The primary sample always gets the alt genotype
                Genotype sourceGt = primaryVc.getGenotype(primarySample)
                Map<String, Object> formatFields = new LinkedHashMap()
                if(sourceGt != null) {
                    Object cr = sourceGt.getExtendedAttribute("CR")
                    Object nc = sourceGt.getExtendedAttribute("NC")
                    if(cr != null) {
                        formatFields.put("CR", cr)
                    }
                    if(nc != null) {
                        formatFields.put("NC", nc)
                    }
                }
                genotypes.add(GenotypeBuilder.create(sample, [refAllele, altAllele], formatFields))
            }
            else {
                // No match - hom-ref
                genotypes.add(GenotypeBuilder.create(sample, [refAllele, refAllele]))
            }
        }

        VariantContextBuilder builder = new VariantContextBuilder()
        builder.chr(primaryVc.getContig())
        builder.start(primaryVc.getStart())
        builder.stop(primaryVc.getEnd())
        builder.log10PError(primaryVc.getLog10PError())
        builder.alleles((Collection<Allele>)primaryVc.getAlleles())
        builder.genotypes((List<Genotype>)genotypes)

        // Copy INFO attributes from primary
        for(Map.Entry<String, Object> attr in primaryVc.getAttributes().entrySet()) {
            builder.attribute(attr.getKey(), attr.getValue())
        }

        // Copy filters from primary
        if(primaryVc.filtersWereApplied()) {
            if(primaryVc.getFilters().isEmpty()) {
                builder.passFilters()
            }
            else {
                builder.filters(primaryVc.getFilters())
            }
        }

        return builder.make()
    }

    /**
     * Write the merged VCF output.
     *
     * @param w output writer
     * @param sourceHeader a VCF header from one of the input files (used for meta-data lines)
     * @param allSamples ordered list of all sample names
     * @param variants list of merged variants to write
     */
    void writeVCF(Writer w, VCFHeader sourceHeader, List<String> allSamples, List<VariantContext> variants) {

        // Build header with all samples
        Set<VCFHeaderLine> metaData = sourceHeader.getMetaDataInInputOrder()
        VCFHeader outputHeader = new VCFHeader(metaData, allSamples)

        // Write header
        w.write(outputHeader.getMetaDataInInputOrder().collect { '##' + it }.join('\n') + '\n')
        w.write((['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + allSamples).join('\t') + '\n')

        // Sort by genomic position
        List<VariantContext> sorted = variants.sort { VariantContext vc ->
            XPos.computePos(vc.getContig(), vc.getStart())
        }

        VCFEncoder encoder = new VCFEncoder(outputHeader, true, true)
        for(VariantContext vc in sorted) {
            encoder.write(w, vc)
            w.write('\n')
        }

        log.info "Wrote ${sorted.size()} variants to output VCF"
    }
}
