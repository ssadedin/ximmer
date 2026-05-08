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
