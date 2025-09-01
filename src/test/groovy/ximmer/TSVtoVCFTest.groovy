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
            columns: ["copy_number", "coverage_ratio"],
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
        
        // Change the implementation and this test so that SVTYPE comes out as DEL in this scenario AI!
        assert variant.getAttribute("SVTYPE") == "DUP,DEL"
        assert variant.alleles.size() == 2
        assert variant.alleles[0].displayString == "G"
        assert variant.alleles[1].displayString == "<DEL>"
    }

}
