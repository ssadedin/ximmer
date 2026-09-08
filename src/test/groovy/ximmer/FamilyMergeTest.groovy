package ximmer

import gngs.Region
import gngs.Regions
import htsjdk.variant.variantcontext.*
import htsjdk.variant.vcf.*
import org.junit.Test

class FamilyMergeTest {

    // ---------------------------------------------------------------
    // Helper to build a minimal VariantContext for testing
    // ---------------------------------------------------------------

    /**
     * Build a test VariantContext with the given properties.
     * Mimics the structure produced by TSVtoVCF.
     */
    private VariantContext makeVariant(Map params) {
        String chr = params.chr ?: 'chr1'
        int start = params.start ?: 1000
        int end = params.end ?: 2000
        String svtype = params.svtype ?: 'DEL'
        String sample = params.sample ?: 'PROBAND'
        List<String> allSamples = params.allSamples ?: [sample]
        double qual = params.qual ?: 30.0d
        int callers = params.callers ?: 1
        double cr = params.cr ?: 0.5d
        int nc = params.nc ?: 1

        String refBase = params.refBase ?: 'A'
        Allele refAllele = Allele.create(refBase, true)
        Allele altAllele = Allele.create('<' + svtype + '>')

        int svLen = (end - start) * (svtype == 'DEL' ? -1 : 1)

        List<Genotype> gts = allSamples.collect { String s ->
            if(s == sample) {
                Map formatFields = [CR: cr, NC: nc]
                return GenotypeBuilder.create(s, [refAllele, altAllele], formatFields)
            }
            return GenotypeBuilder.create(s, [refAllele, refAllele])
        }

        return new VariantContextBuilder()
                .chr(chr)
                .start(start)
                .stop(end)
                .log10PError(-qual / 10.0d)
                .attribute("SVTYPE", svtype)
                .attribute("END", end)
                .attribute("SVLEN", svLen)
                .attribute("CR", cr)
                .attribute("CN", svtype == 'DEL' ? 1 : 3)
                .attribute("CALLERS", callers)
                .attribute("CALLEDBY", ['testcaller'])
                .alleles([refAllele, altAllele])
                .genotypes(gts)
                .passFilters()
                .make()
    }

    // ---------------------------------------------------------------
    // Tests for exactMatchKey
    // ---------------------------------------------------------------

    @Test
    void testExactMatchKeySameVariants() {
        def vc1 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DEL', sample: 'S1')
        def vc2 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DEL', sample: 'S2')

        assert FamilyMerge.exactMatchKey(vc1) == FamilyMerge.exactMatchKey(vc2)
    }

    @Test
    void testExactMatchKeyDifferentPosition() {
        def vc1 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DEL')
        def vc2 = makeVariant(chr: 'chr1', start: 1000, end: 3000, svtype: 'DEL')

        assert FamilyMerge.exactMatchKey(vc1) != FamilyMerge.exactMatchKey(vc2)
    }

    @Test
    void testExactMatchKeyDifferentType() {
        def vc1 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DEL')
        def vc2 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DUP')

        assert FamilyMerge.exactMatchKey(vc1) != FamilyMerge.exactMatchKey(vc2)
    }

    @Test
    void testExactMatchKeyDifferentChrom() {
        def vc1 = makeVariant(chr: 'chr1', start: 1000, end: 2000, svtype: 'DEL')
        def vc2 = makeVariant(chr: 'chr2', start: 1000, end: 2000, svtype: 'DEL')

        assert FamilyMerge.exactMatchKey(vc1) != FamilyMerge.exactMatchKey(vc2)
    }

    // ---------------------------------------------------------------
    // Tests for findNonRefSample
    // ---------------------------------------------------------------

    @Test
    void testFindNonRefSample() {
        def vc = makeVariant(sample: 'CHILD', allSamples: ['CHILD', 'MOM', 'DAD'])
        assert FamilyMerge.findNonRefSample(vc) == 'CHILD'
    }

    @Test
    void testFindNonRefSampleAllRef() {
        Allele ref = Allele.create('A', true)
        List<Genotype> gts = ['S1', 'S2'].collect {
            GenotypeBuilder.create(it, [ref, ref])
        }
        def vc = new VariantContextBuilder()
                .chr('chr1').start(100).stop(200)
                .alleles([ref, Allele.create('<DEL>')])
                .genotypes(gts)
                .make()

        assert FamilyMerge.findNonRefSample(vc) == null
    }

    // ---------------------------------------------------------------
    // Tests for buildJointVariant
    // ---------------------------------------------------------------

    @Test
    void testBuildJointVariantProbandOnly() {
        def merger = new FamilyMerge()
        def allSamples = ['PROBAND', 'MOM', 'DAD']

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'], cr: 0.3d, nc: 2)

        // No parent matches
        Map<String, VariantContext> matches = [:]
        def joint = merger.buildJointVariant(probandVc, 'PROBAND', allSamples, matches)

        // Proband should be het
        def probandGt = joint.getGenotype('PROBAND')
        assert probandGt != null
        assert !probandGt.isHomRef()
        assert probandGt.getExtendedAttribute("CR") == 0.3d
        assert probandGt.getExtendedAttribute("NC") == 2

        // Parents should be hom-ref
        assert joint.getGenotype('MOM').isHomRef()
        assert joint.getGenotype('DAD').isHomRef()

        // Coordinates from primary
        assert joint.contig == 'chr1'
        assert joint.start == 1000
        assert joint.end == 2000
        assert joint.getAttribute("SVTYPE") == 'DEL'
    }

    @Test
    void testBuildJointVariantWithParentMatch() {
        def merger = new FamilyMerge()
        def allSamples = ['PROBAND', 'MOM', 'DAD']

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'], cr: 0.2d, nc: 3)
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'], cr: 0.4d, nc: 1)

        Map<String, VariantContext> matches = [MOM: momVc]
        def joint = merger.buildJointVariant(probandVc, 'PROBAND', allSamples, matches)

        // Proband should be het with own CR/NC
        def probandGt = joint.getGenotype('PROBAND')
        assert !probandGt.isHomRef()
        assert probandGt.getExtendedAttribute("CR") == 0.2d
        assert probandGt.getExtendedAttribute("NC") == 3

        // MOM should be het with her own CR/NC
        def momGt = joint.getGenotype('MOM')
        assert !momGt.isHomRef()
        assert momGt.getExtendedAttribute("CR") == 0.4d
        assert momGt.getExtendedAttribute("NC") == 1

        // DAD should be hom-ref
        assert joint.getGenotype('DAD').isHomRef()
    }

    @Test
    void testBuildJointVariantBothParentsMatch() {
        def merger = new FamilyMerge()
        def allSamples = ['PROBAND', 'MOM', 'DAD']

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'], cr: 0.1d, nc: 2)
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'], cr: 0.3d, nc: 1)
        def dadVc = makeVariant(sample: 'DAD', allSamples: ['DAD'], cr: 0.5d, nc: 1)

        Map<String, VariantContext> matches = [MOM: momVc, DAD: dadVc]
        def joint = merger.buildJointVariant(probandVc, 'PROBAND', allSamples, matches)

        // All three should be het
        assert !joint.getGenotype('PROBAND').isHomRef()
        assert !joint.getGenotype('MOM').isHomRef()
        assert !joint.getGenotype('DAD').isHomRef()

        // Each has their own CR
        assert joint.getGenotype('PROBAND').getExtendedAttribute("CR") == 0.1d
        assert joint.getGenotype('MOM').getExtendedAttribute("CR") == 0.3d
        assert joint.getGenotype('DAD').getExtendedAttribute("CR") == 0.5d
    }

    // ---------------------------------------------------------------
    // Tests for mergeVariants (full trio merge logic)
    // ---------------------------------------------------------------

    @Test
    void testMergeProbandEventMatchedByOneParent() {
        def merger = new FamilyMerge()

        // Same exact coordinates
        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'])
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'])

        Map<String, List<VariantContext>> variants = [
            PROBAND: [probandVc],
            MOM: [momVc],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // Should produce exactly 1 merged variant
        assert result.size() == 1

        def merged = result[0]
        assert !merged.getGenotype('PROBAND').isHomRef()
        assert !merged.getGenotype('MOM').isHomRef()
        assert merged.getGenotype('DAD').isHomRef()
    }

    @Test
    void testMergeProbandEventNoParentMatch() {
        def merger = new FamilyMerge()

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'],
                                     chr: 'chr1', start: 1000, end: 2000)
        // Parent has event at different position
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 chr: 'chr1', start: 5000, end: 6000)

        Map<String, List<VariantContext>> variants = [
            PROBAND: [probandVc],
            MOM: [momVc],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // 1 proband event + 1 unconsumed parent event = 2
        assert result.size() == 2

        // Find the proband-driven event
        def probandEvent = result.find { !it.getGenotype('PROBAND').isHomRef() }
        assert probandEvent != null
        assert probandEvent.getGenotype('MOM').isHomRef()
        assert probandEvent.getGenotype('DAD').isHomRef()

        // Find the parent-only event
        def parentEvent = result.find { it.getGenotype('PROBAND').isHomRef() }
        assert parentEvent != null
        assert !parentEvent.getGenotype('MOM').isHomRef()
    }

    @Test
    void testMergeParentOnlyEventsMatchEachOther() {
        def merger = new FamilyMerge()

        // No proband events; both parents have same event
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 chr: 'chr2', start: 3000, end: 4000)
        def dadVc = makeVariant(sample: 'DAD', allSamples: ['DAD'],
                                 chr: 'chr2', start: 3000, end: 4000)

        Map<String, List<VariantContext>> variants = [
            PROBAND: [],
            MOM: [momVc],
            DAD: [dadVc]
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // Should produce 1 merged variant (both parents matched)
        assert result.size() == 1

        def merged = result[0]
        assert merged.getGenotype('PROBAND').isHomRef()
        assert !merged.getGenotype('MOM').isHomRef()
        assert !merged.getGenotype('DAD').isHomRef()
    }

    @Test
    void testMergeParentOnlyEventsDifferentPositions() {
        def merger = new FamilyMerge()

        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 chr: 'chr1', start: 1000, end: 2000)
        def dadVc = makeVariant(sample: 'DAD', allSamples: ['DAD'],
                                 chr: 'chr1', start: 5000, end: 6000)

        Map<String, List<VariantContext>> variants = [
            PROBAND: [],
            MOM: [momVc],
            DAD: [dadVc]
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // 2 separate parent-only events
        assert result.size() == 2

        // Both should have proband as hom-ref
        assert result.every { it.getGenotype('PROBAND').isHomRef() }
    }

    @Test
    void testMergeSvtypeMismatchPreventsMatch() {
        def merger = new FamilyMerge()

        // Same position but different SVTYPE
        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'],
                                     svtype: 'DEL')
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 svtype: 'DUP')

        Map<String, List<VariantContext>> variants = [
            PROBAND: [probandVc],
            MOM: [momVc],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // 2 separate events (no match due to type mismatch)
        assert result.size() == 2

        // Proband event should have MOM as hom-ref
        def probandEvent = result.find { !it.getGenotype('PROBAND').isHomRef() }
        assert probandEvent.getGenotype('MOM').isHomRef()
    }

    @Test
    void testMergeMultipleEventsOnSameChromosome() {
        def merger = new FamilyMerge()

        // Proband has 2 events, parent matches one
        def probandVc1 = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'],
                                      chr: 'chr1', start: 1000, end: 2000)
        def probandVc2 = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'],
                                      chr: 'chr1', start: 5000, end: 6000)
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 chr: 'chr1', start: 1000, end: 2000)

        Map<String, List<VariantContext>> variants = [
            PROBAND: [probandVc1, probandVc2],
            MOM: [momVc],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // 2 proband events, no leftover parent events
        assert result.size() == 2

        // The event at 1000-2000 should have MOM as het
        def matchedEvent = result.find { it.start == 1000 }
        assert !matchedEvent.getGenotype('MOM').isHomRef()

        // The event at 5000-6000 should have MOM as hom-ref
        def unmatchedEvent = result.find { it.start == 5000 }
        assert unmatchedEvent.getGenotype('MOM').isHomRef()
    }

    @Test
    void testMergeSingleParent() {
        def merger = new FamilyMerge()

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'])
        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'])

        Map<String, List<VariantContext>> variants = [
            PROBAND: [probandVc],
            MOM: [momVc]
        ]

        // Only one parent
        def result = merger.mergeVariants('PROBAND', ['MOM'], variants)

        assert result.size() == 1
        assert !result[0].getGenotype('PROBAND').isHomRef()
        assert !result[0].getGenotype('MOM').isHomRef()
    }

    @Test
    void testMergePreservesInfoAttributes() {
        def merger = new FamilyMerge()
        def allSamples = ['PROBAND', 'MOM']

        def probandVc = makeVariant(sample: 'PROBAND', allSamples: ['PROBAND'],
                                     chr: 'chr3', start: 2000, end: 3000, svtype: 'DEL',
                                     callers: 3, cr: 0.15d)

        Map<String, VariantContext> matches = [:]
        def joint = merger.buildJointVariant(probandVc, 'PROBAND', allSamples, matches)

        // INFO fields should be carried from primary
        assert joint.getAttribute("SVTYPE") == 'DEL'
        assert joint.getAttribute("END") == 3000
        assert joint.getAttribute("SVLEN") == -1000
        assert joint.getAttribute("CALLERS") == 3
        assert joint.getAttribute("CR") == 0.15d
    }

    @Test
    void testMergePreservesFilters() {
        def merger = new FamilyMerge()
        def allSamples = ['PROBAND', 'MOM']

        // Build a variant with a filter applied
        Allele ref = Allele.create('A', true)
        Allele alt = Allele.create('<DEL>')
        def filteredVc = new VariantContextBuilder()
                .chr('chr1').start(1000).stop(2000)
                .attribute("SVTYPE", "DEL")
                .attribute("END", 2000)
                .attribute("SVLEN", -1000)
                .attribute("CR", 0.5d)
                .attribute("CN", 1)
                .attribute("CALLERS", 1)
                .attribute("CALLEDBY", ['xhmm'])
                .alleles([ref, alt])
                .genotypes([GenotypeBuilder.create('PROBAND', [ref, alt], [CR: 0.5d, NC: 1])])
                .filters(['LOW_CALLERS'] as Set)
                .make()

        Map<String, VariantContext> matches = [:]
        def joint = merger.buildJointVariant(filteredVc, 'PROBAND', allSamples, matches)

        assert joint.getFilters().contains('LOW_CALLERS')
    }

    @Test
    void testMergeNoVariantsForProband() {
        def merger = new FamilyMerge()

        def momVc = makeVariant(sample: 'MOM', allSamples: ['MOM'],
                                 chr: 'chr1', start: 1000, end: 2000)

        Map<String, List<VariantContext>> variants = [
            PROBAND: [],
            MOM: [momVc],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        // 1 parent-only event
        assert result.size() == 1
        assert result[0].getGenotype('PROBAND').isHomRef()
        assert !result[0].getGenotype('MOM').isHomRef()
        assert result[0].getGenotype('DAD').isHomRef()
    }

    @Test
    void testMergeEmptyInput() {
        def merger = new FamilyMerge()

        Map<String, List<VariantContext>> variants = [
            PROBAND: [],
            MOM: [],
            DAD: []
        ]

        def result = merger.mergeVariants('PROBAND', ['MOM', 'DAD'], variants)

        assert result.isEmpty()
    }
}
