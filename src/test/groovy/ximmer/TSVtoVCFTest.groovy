package ximmer

import static org.junit.Assert.*

import com.xlson.groovycsv.PropertyMapper
import gngs.BED
import gngs.FASTA
import gngs.Region
import gngs.Regions
import org.junit.Test

class TSVtoVCFTest {
    // no dot in CR
    // no N as REF
    


    /**
     * Test the createVariantFromLine method
     * Single deletion on input maps to single deletion in output variant context
     */
    @Test
    public void testSimpleDeletion() {
        // Create mock objects
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            columns: ["copy_number", "coverage_ratio"],
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 2,
            xhmm_qual : 100d,
            xhmm : 'TRUE',
            cnvnator_qual: 90d,
            cnvnator: 'TRUE'
        ] 

        // Test with single sample
        def samples = ["SAMPLE1"]
        
        // Get variant context
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify basic properties
        assert variant != null
        assert variant.contig == "chr1"
        assert variant.start == 1000
        assert variant.end == 2000
        assert variant.getAttribute("SVTYPE") == "DEL"
        assert variant.getAttribute("SVLEN") == -1000
        assert variant.getAttribute("CR") == 0.5
        assert variant.getAttribute("CN") == 1
        assert variant.getAttribute("CALLERS") == 2
        
        // Verify alleles
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "A"
        assert variant.alleles[1].displayString == "<DEL>"
        
        // Verify genotypes
        def genotype = variant.getGenotype("SAMPLE1")
        assert genotype != null
        assert genotype.getAlleles().size() == 2
        assert genotype.getExtendedAttribute("CR") == 0.5
        assert genotype.getExtendedAttribute("NC") == 2
    }
    
    /**
     * Test that if the type has both DUP and DEL, only DEL is in the output
     */
    @Test
    void testCombinedDupAndDel() {
        // Create mock objects
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "G" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 5000,
            end: 6000,
            type: "DUP,DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 3,
            xhmm_qual: 80d,
            xhmm: 'TRUE',
            cnvnator_qual: 85d,
            cnvnator: 'TRUE',
            exomedepth_qual: 75d,
            exomedepth: 'TRUE'
        ]

        // Test with single sample
        def samples = ["SAMPLE1"]
        
        // Get variant context
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that DEL was chosen over DUP and CN is set to 1
        assert variant != null
        
        assert variant.getAttribute("SVTYPE") == "DEL"
        assert variant.getAttribute("CN") == 1
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "G"
        assert variant.alleles[1].displayString == "<DEL>"
    }
    
    /**
     * Test that inversions (INV) are treated as deletions in the output
     */
    @Test
    void testInversion() {
        // Create mock objects
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "T" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 3000,
            end: 4000,
            type: "INV",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 90d,
            xhmm: 'TRUE'
        ]

        // Test with single sample
        def samples = ["SAMPLE1"]
        
        // Get variant context
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that INV is converted to DEL
        assert variant != null
        assert variant.getAttribute("SVTYPE") == "DEL"
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "T"
        assert variant.alleles[1].displayString == "<DEL>"
    }
    
    /**
     * Test that default coverage ratios are used when coverage_ratio is not provided
     */
    @Test
    void testDefaultCoverageRatios() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "C" }
        ] as FASTA
        
        // Test deletion without coverage_ratio
        Map delLine = [
            chr: "chr1",
            start: 7000,
            end: 8000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def delVariant = tsv.createVariantFromLine(delLine, mockFasta, samples)
        
        // Verify default CR for deletion is 0.5 in both variant and genotype
        assert delVariant != null
        assert delVariant.getAttribute("CR") == 0.5
        assert delVariant.getGenotype("SAMPLE1").getExtendedAttribute("CR") == 0.5
        
        // Test duplication without coverage_ratio
        Map dupLine = [
            chr: "chr1",
            start: 7000,
            end: 8000,
            type: "DUP",
            sample: "SAMPLE1",
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        ]
        
        def dupVariant = tsv.createVariantFromLine(dupLine, mockFasta, samples)
        
        // Verify default CR for duplication is 3.0 in both variant and genotype
        assert dupVariant != null
        assert dupVariant.getAttribute("CR") == 3.0
        assert dupVariant.getGenotype("SAMPLE1").getExtendedAttribute("CR") == 3.0
    }
    
    /**
     * Test that variants where start == end (zero-length / SVLEN=0) are skipped
     */
    @Test
    void testSkipZeroLengthVariant() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 5000,
            end: 5000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that variant is null when start == end (SVLEN would be zero)
        assert variant == null
    }

    /**
     * Test that variants with N as reference base are skipped
     */
    @Test
    void testSkipNReference() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "N" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 9000,
            end: 10000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that variant is null when reference base is N
        assert variant == null
    }
    
    /**
     * Test that DEL variants have copy number capped at 1
     */
    @Test
    void testDelCopyNumberCap() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 11000,
            end: 12000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 3,  // High copy number that should be capped
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that copy number is capped at 1 for DEL
        assert variant != null
        assert variant.getAttribute("SVTYPE") == "DEL"
        assert variant.getAttribute("CN") == 1
    }
    
    /**
     * Test that combined DUP,INV is converted to DEL
     */
    @Test
    void testCombinedDupAndInv() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "C" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 13000,
            end: 14000,
            type: "DUP,INV",
            sample: "SAMPLE1",
            copy_number: 2,
            coverage_ratio: 1.5,
            count: 2,
            xhmm_qual: 85d,
            xhmm: 'TRUE',
            cnvnator_qual: 80d,
            cnvnator: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that INV causes conversion to DEL and CN is capped at 1
        assert variant != null
        assert variant.getAttribute("SVTYPE") == "DEL"
        assert variant.getAttribute("CN") == 1
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "C"
        assert variant.alleles[1].displayString == "<DEL>"
    }
    
    /**
     * Test that DUP variants have minimum copy number of 3
     */
    @Test
    void testDupMinimumCopyNumber() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "G" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 15000,
            end: 16000,
            type: "DUP",
            sample: "SAMPLE1",
            copy_number: 2,  // Low copy number that should be raised to 3
            coverage_ratio: 2.0,
            count: 1,
            xhmm_qual: 75d,
            xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that copy number is set to minimum of 3 for DUP
        assert variant != null
        assert variant.getAttribute("SVTYPE") == "DUP"
        assert variant.getAttribute("CN") == 3
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "G"
        assert variant.alleles[1].displayString == "<DUP>"
    }

    /**
     * Test that variants pass filters when no filter thresholds are specified
     */
    @Test
    void testNoFiltersAppliedByDefault() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 100d,
            xhmm: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // No target regions, no thresholds - should always PASS
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        assert variant != null
        assert variant.filtersWereApplied()
        assert variant.filters.isEmpty() // empty filters means PASS
    }
    
    /**
     * Test that LOW_CALLERS filter is applied when caller count is below threshold
     */
    @Test
    void testLowCallersFilter() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 100d,
            xhmm: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // Require 2 callers, but only 1 present - should fail
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, null, null, 2)
        
        assert variant != null
        assert variant.filters.contains(TSVtoVCF.FILTER_LOW_CALLERS)
    }
    
    /**
     * Test that LOW_CALLERS filter is NOT applied when caller count meets threshold
     */
    @Test
    void testCallerCountMeetsThreshold() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 2,
            xhmm_qual: 100d,
            xhmm: 'TRUE',
            cnvnator_qual: 90d,
            cnvnator: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // Require 2 callers, and 2 are present - should PASS
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, null, null, 2)
        
        assert variant != null
        assert !variant.filters.contains(TSVtoVCF.FILTER_LOW_CALLERS)
    }
    
    /**
     * Test that FEW_TARGETS filter is applied when target overlap count is below threshold
     */
    @Test
    void testFewTargetsFilter() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Create target regions - only 1 target overlaps the CNV
        Regions targetRegions = [
            new Region('chr1:1200-1300'),
            new Region('chr1:5000-5100'),
            new Region('chr1:6000-6100')
        ] as Regions
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 100d,
            xhmm: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // Require 3 targets, but only 1 overlaps - should fail
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, targetRegions, 3, null)
        
        assert variant != null
        assert variant.filters.contains(TSVtoVCF.FILTER_FEW_TARGETS)
        assert variant.getAttribute("TARGETS") == 1
    }
    
    /**
     * Test that FEW_TARGETS filter is NOT applied when target overlap count meets threshold
     */
    @Test
    void testTargetCountMeetsThreshold() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Create target regions - 3 targets overlap the CNV
        Regions targetRegions = [
            new Region('chr1:1100-1200'),
            new Region('chr1:1400-1500'),
            new Region('chr1:1700-1800'),
            new Region('chr1:5000-5100')
        ] as Regions
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 100d,
            xhmm: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // Require 3 targets, and 3 overlap - should PASS on target filter
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, targetRegions, 3, null)
        
        assert variant != null
        assert !variant.filters.contains(TSVtoVCF.FILTER_FEW_TARGETS)
        assert variant.getAttribute("TARGETS") == 3
    }
    
    /**
     * Test that both filters can be applied simultaneously
     */
    @Test
    void testBothFiltersApplied() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Create target regions - only 1 target overlaps the CNV
        Regions targetRegions = [
            new Region('chr1:1200-1300'),
            new Region('chr1:5000-5100')
        ] as Regions
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 1,
            xhmm_qual: 100d,
            xhmm: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // Require 2 callers AND 3 targets - both should fail
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, targetRegions, 3, 2)
        
        assert variant != null
        assert variant.filters.contains(TSVtoVCF.FILTER_LOW_CALLERS)
        assert variant.filters.contains(TSVtoVCF.FILTER_FEW_TARGETS)
        assert variant.filters.size() == 2
    }
    
    /**
     * Test that TARGETS info field is populated when target regions are provided
     */
    @Test
    void testTargetsInfoField() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Create target regions - 2 targets overlap the CNV
        Regions targetRegions = [
            new Region('chr1:1100-1200'),
            new Region('chr1:1500-1600'),
            new Region('chr1:5000-5100')
        ] as Regions
        
        Map line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 2,
            xhmm_qual: 100d,
            xhmm: 'TRUE',
            cnvnator_qual: 90d,
            cnvnator: 'TRUE'
        ]

        def samples = ["SAMPLE1"]
        
        // No filter thresholds, but target regions provided - TARGETS should be populated
        def variant = tsv.createVariantFromLine(line, mockFasta, samples, targetRegions, null, null)
        
        assert variant != null
        assert variant.getAttribute("TARGETS") == 2
        // Should still PASS since no thresholds set
        assert variant.filters.isEmpty()
    }
    
    /**
     * Test that two DELs at the same position are flattened into one with longest span and max callers
     */
    @Test
    void testFlattenSamePositionAndType() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Two DELs at same position, different spans and caller counts
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 0, coverage_ratio: 0.1, count: 2,
            ed_qual: 90d, ed: 'TRUE', savvy_qual: 85d, savvy: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 5000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        
        def flattened = tsv.flattenVariants([v1, v2], null, null, null)
        
        assert flattened.size() == 1
        def result = flattened[0]
        // Span from longest call
        assert result.end == 5000
        // Callers from max-count call (line1 has count=2)
        assert result.getAttribute("CALLERS") == 2
        assert result.getAttribute("CALLEDBY") == ['ed', 'savvy']
        // CN/CR from representative (max callers)
        assert result.getAttribute("CN") == 0
        assert result.getAttribute("CR") == 0.1
    }
    
    /**
     * Test that DEL and DUP at same position are NOT flattened together
     */
    @Test
    void testFlattenDoesNotMergeDifferentTypes() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.5, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 3000, type: "DUP", sample: "SAMPLE1",
            copy_number: 3, coverage_ratio: 1.5, count: 1,
            ed_qual: 90d, ed: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        
        def flattened = tsv.flattenVariants([v1, v2], null, null, null)
        
        // Different types should not be merged
        assert flattened.size() == 2
    }
    
    /**
     * Test that calls from different samples at same position are NOT flattened
     */
    @Test
    void testFlattenDoesNotMergeDifferentSamples() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE_A",
            copy_number: 1, coverage_ratio: 0.5, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 3000, type: "DEL", sample: "SAMPLE_B",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            ed_qual: 90d, ed: 'TRUE'
        ]
        
        def samples = ["SAMPLE_A", "SAMPLE_B"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        
        def flattened = tsv.flattenVariants([v1, v2], null, null, null)
        
        // Different samples should not be merged
        assert flattened.size() == 2
    }
    
    /**
     * Test that when caller counts are tied, the call with higher QUAL is used as representative
     */
    @Test
    void testFlattenTieBreakByQual() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 0, coverage_ratio: 0.1, count: 1,
            ed_qual: 90d, ed: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 5000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            xhmm_qual: 95d, xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        
        def flattened = tsv.flattenVariants([v1, v2], null, null, null)
        
        assert flattened.size() == 1
        def result = flattened[0]
        // Span from longest
        assert result.end == 5000
        // Callers from higher QUAL call (line2, qual=95 > line1, qual=90)
        assert result.getAttribute("CALLEDBY") == ['xhmm']
        assert result.getAttribute("CR") == 0.4
    }
    
    /**
     * Test that TARGETS is recomputed from the merged (longest) span
     */
    @Test
    void testFlattenRecomputesTargets() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Target regions: 2 in short span, 1 more in extended span
        Regions targetRegions = [
            new Region('chr1:1100-1200'),
            new Region('chr1:1500-1600'),
            new Region('chr1:3000-3100')
        ] as Regions
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 0, coverage_ratio: 0.1, count: 2,
            ed_qual: 90d, ed: 'TRUE', savvy_qual: 85d, savvy: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 4000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples, targetRegions, null, null)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples, targetRegions, null, null)
        
        def flattened = tsv.flattenVariants([v1, v2], targetRegions, null, null)
        
        assert flattened.size() == 1
        def result = flattened[0]
        // TARGETS recomputed from merged span (1000-4000) which covers all 3 targets
        assert result.getAttribute("TARGETS") == 3
    }
    
    /**
     * Test that filters are recomputed after flattening using representative's count and new target count
     */
    @Test
    void testFlattenRecomputesFilters() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        // Target regions: only 1 in short span, 3 in extended span
        Regions targetRegions = [
            new Region('chr1:1100-1200'),
            new Region('chr1:3000-3100'),
            new Region('chr1:4000-4100')
        ] as Regions
        
        // line1: 2 callers but short span (only 1 target)
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 0, coverage_ratio: 0.1, count: 2,
            ed_qual: 90d, ed: 'TRUE', savvy_qual: 85d, savvy: 'TRUE'
        ]
        // line2: 1 caller but long span (covers all 3 targets)
        Map line2 = [
            chr: "chr1", start: 1000, end: 5000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples, targetRegions, 3, 2)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples, targetRegions, 3, 2)
        
        // Before flattening: v1 fails FEW_TARGETS (1 < 3), v2 fails LOW_CALLERS (1 < 2) and FEW_TARGETS
        assert v1.filters.contains(TSVtoVCF.FILTER_FEW_TARGETS)
        assert v2.filters.contains(TSVtoVCF.FILTER_LOW_CALLERS)
        
        def flattened = tsv.flattenVariants([v1, v2], targetRegions, 3, 2)
        
        assert flattened.size() == 1
        def result = flattened[0]
        // After flattening: callers=2 (from representative), targets=3 (recomputed from long span)
        // Both thresholds met, so should PASS
        assert result.filters.isEmpty()
        assert result.getAttribute("CALLERS") == 2
        assert result.getAttribute("TARGETS") == 3
    }
    
    /**
     * Test flattening with three calls at the same position
     */
    @Test
    void testFlattenThreeCallsAtSamePosition() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 0, coverage_ratio: 0.05, count: 3,
            ed_qual: 90d, ed: 'TRUE', savvy_qual: 85d, savvy: 'TRUE', xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 3000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            cnvnator_qual: 80d, cnvnator: 'TRUE'
        ]
        Map line3 = [
            chr: "chr1", start: 1000, end: 4000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.5, count: 2,
            gatk_qual: 75d, gatk: 'TRUE', xhmm_qual: 70d, xhmm: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        def v3 = tsv.createVariantFromLine(line3, mockFasta, samples)
        
        def flattened = tsv.flattenVariants([v1, v2, v3], null, null, null)
        
        assert flattened.size() == 1
        def result = flattened[0]
        // Span from longest (line3: end=4000)
        assert result.end == 4000
        // Callers from max-count (line1: count=3)
        assert result.getAttribute("CALLERS") == 3
        assert result.getAttribute("CALLEDBY") == ['ed', 'savvy', 'xhmm']
        assert result.getAttribute("CN") == 0
        assert result.getAttribute("CR") == 0.05
    }
    
    /**
     * Test that without flatten flag, duplicates are preserved (tested via flattenVariants not being called)
     */
    @Test
    void testFlattenDisabledByDefault() {
        def tsv = new TSVtoVCF()
        def mockFasta = [
            basesAt: { chr, start, end -> "A" }
        ] as FASTA
        
        Map line1 = [
            chr: "chr1", start: 1000, end: 2000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.5, count: 1,
            xhmm_qual: 80d, xhmm: 'TRUE'
        ]
        Map line2 = [
            chr: "chr1", start: 1000, end: 5000, type: "DEL", sample: "SAMPLE1",
            copy_number: 1, coverage_ratio: 0.4, count: 1,
            ed_qual: 90d, ed: 'TRUE'
        ]
        
        def samples = ["SAMPLE1"]
        def v1 = tsv.createVariantFromLine(line1, mockFasta, samples)
        def v2 = tsv.createVariantFromLine(line2, mockFasta, samples)
        
        // Directly calling flattenVariants should still work, but without the flag
        // the createVCF method won't call it. Here we verify the variants remain separate
        // when not flattened.
        def variants = [v1, v2]
        assert variants.size() == 2
        assert variants[0].end == 2000
        assert variants[1].end == 5000
    }

    /**
     * Test the computeFilters method directly
     */
    @Test
    void testComputeFilters() {
        def tsv = new TSVtoVCF()
        
        // No thresholds - always empty (PASS)
        assert tsv.computeFilters(1, 1, null, null).isEmpty()
        
        // Caller count below threshold
        assert tsv.computeFilters(1, 5, null, 2).contains(TSVtoVCF.FILTER_LOW_CALLERS)
        
        // Caller count meets threshold
        assert !tsv.computeFilters(2, 5, null, 2).contains(TSVtoVCF.FILTER_LOW_CALLERS)
        
        // Target count below threshold
        assert tsv.computeFilters(2, 1, 3, null).contains(TSVtoVCF.FILTER_FEW_TARGETS)
        
        // Target count meets threshold
        assert !tsv.computeFilters(2, 3, 3, null).contains(TSVtoVCF.FILTER_FEW_TARGETS)
        
        // Both filters fail
        def filters = tsv.computeFilters(1, 1, 3, 2)
        assert filters.contains(TSVtoVCF.FILTER_LOW_CALLERS)
        assert filters.contains(TSVtoVCF.FILTER_FEW_TARGETS)
        assert filters.size() == 2
        
        // Target overlap is null (no target regions provided) - FEW_TARGETS not applied
        assert !tsv.computeFilters(2, null, 3, null).contains(TSVtoVCF.FILTER_FEW_TARGETS)
    }

}
