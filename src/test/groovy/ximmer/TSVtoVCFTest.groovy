package ximmer

import static org.junit.Assert.*

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
        
        // Create mock input line
        def line = [
            chr: "chr1",
            start: 1000,
            end: 2000,
            type: "DEL",
            sample: "SAMPLE1",
            columns: ["copy_number", "coverage_ratio"],
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 2
        ] as PropertyMapper
        
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

}
