package ximmer

import static org.junit.Assert.*

import com.xlson.groovycsv.PropertyMapper
import gngs.FASTA
import org.junit.Test

class TSVtoVCFTest {
    // no dot in CR
    // no N as REF
    
    /**
     * Creates a PropertyMapper from a map of values
     */
    private PropertyMapper createPropertyMapper(Map lineValues) {
        PropertyMapper mapper = new PropertyMapper()
        mapper.columns = lineValues*.key.indexed().collectEntries { [it.value, it.key] }
        mapper.values = lineValues*.value
        return mapper
    }

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
        
        PropertyMapper line = createPropertyMapper(
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
            xhmm : 'TRUE'
        ) 

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
        
        PropertyMapper line = createPropertyMapper(
            chr: "chr1",
            start: 5000,
            end: 6000,
            type: "DUP,DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            coverage_ratio: 0.5,
            count: 3,
            xhmm_qual: 80d,
            xhmm: 'TRUE'
        )

        // Test with single sample
        def samples = ["SAMPLE1"]
        
        // Get variant context
        def variant = tsv.createVariantFromLine(line, mockFasta, samples)
        
        // Verify that DEL was chosen over DUP
        assert variant != null
        
        assert variant.getAttribute("SVTYPE") == "DEL"
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
        
        PropertyMapper line = createPropertyMapper(
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
        )

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
        PropertyMapper delLine = createPropertyMapper(
            chr: "chr1",
            start: 7000,
            end: 8000,
            type: "DEL",
            sample: "SAMPLE1",
            copy_number: 1,
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        )
        
        def samples = ["SAMPLE1"]
        def delVariant = tsv.createVariantFromLine(delLine, mockFasta, samples)
        
        // Verify default CR for deletion is 0.5 in both variant and genotype
        assert delVariant != null
        assert delVariant.getAttribute("CR") == 0.5
        assert delVariant.getGenotype("SAMPLE1").getExtendedAttribute("CR") == 0.5
        
        // Test duplication without coverage_ratio
        PropertyMapper dupLine = createPropertyMapper(
            chr: "chr1",
            start: 7000,
            end: 8000,
            type: "DUP",
            sample: "SAMPLE1",
            count: 1,
            xhmm_qual: 70d,
            xhmm: 'TRUE'
        )
        
        def dupVariant = tsv.createVariantFromLine(dupLine, mockFasta, samples)
        
        // Verify default CR for duplication is 3.0 in both variant and genotype
        assert dupVariant != null
        assert dupVariant.getAttribute("CR") == 3.0
        assert dupVariant.getGenotype("SAMPLE1").getExtendedAttribute("CR") == 3.0
    }

}
