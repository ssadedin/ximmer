package ximmer.results

import static org.junit.Assert.*

import org.junit.Test
import groovy.json.JsonSlurper

class CNVNatorResultsTest {

	static String cnvnatorResultFile = 'src/test/data/cnvnator/test.cnvnator.vcf'

	@Test
	public void 'TestLoadFileSuccess'() {
		CNVNatorResults c = new CNVNatorResults(cnvnatorResultFile)
	}
	
	@Test(expected = java.io.FileNotFoundException.class)
	public void 'TestLoadFileFail'() {
		CNVNatorResults c = new CNVNatorResults(cnvnatorResultFile + 'blah')
	}
	
	@Test
	public void 'Set sample explicitly'() {
		CNVNatorResults c = new CNVNatorResults(cnvnatorResultFile, 'someSample')
		assert c.sample == 'someSample'
	}

	@Test
	public void "toJson Test"() {
		CNVNatorResults c = new CNVNatorResults(cnvnatorResultFile, 'someSample')
		def b = c.toListMap()
		def jsonStr = """[[chr:chrM, start:1, end:16601, type:DUP], [chr:chr1, start:1, end:10001, type:DEL], [chr:chr1, start:10901, end:26601, type:DUP], [chr:chr1, start:16201, end:37701, type:DUP], [chr:chr1, start:66001, end:133701, type:DEL], [chr:chr1, start:68701, end:162101, type:DEL], [chr:chr1, start:130501, end:265101, type:DUP], [chr:chr1, start:137301, end:276601, type:DUP], [chr:chr1, start:145501, end:298601, type:DEL], [chr:chr1, start:153601, end:309701, type:DEL], [chr:chr1, start:177101, end:404801, type:DEL], [chr:chr1, start:228001, end:467701, type:DUP], [chr:chr1, start:267301, end:585701, type:DEL], [chr:chr1, start:351501, end:737501, type:DEL], [chr:chr1, start:426301, end:875201, type:DUP], [chr:chr1, start:454201, end:918101, type:DUP], [chr:chr1, start:465301, end:936601, type:DUP], [chr:chr1, start:471301, end:992801, type:DEL], [chr:chr1, start:528201, end:1114201, type:DUP], [chr:chr1, start:657901, end:1320401, type:DUP], [chr:chr1, start:672101, end:1348401, type:DUP], [chr:chr1, start:702601, end:1420501, type:DUP], [chr:chr1, start:754301, end:1523401, type:DUP], [chr:chr1, start:793401, end:1637201, type:DUP], [chr:chr1, start:845101, end:1714501, type:DUP], [chr:chr1, start:870301, end:1782701, type:DUP], [chr:chr1, start:913401, end:1925101, type:DUP], [chr:chr1, start:1012801, end:2087501, type:DUP]]
"""
		def jsonSlurper = new JsonSlurper()
		def a = jsonSlurper.parseText(jsonStr)
		//assert a.size() == b.size()
		//assert b[0].sample == a[0].sample
		//assert b[-1].sample == a[-1].sample
	}

}
