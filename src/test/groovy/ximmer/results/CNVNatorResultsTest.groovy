package ximmer.results

import static org.junit.Assert.*

import org.junit.Test
import groovy.json.JsonSlurper

class CNVNatorResultsTest {

	static String cnvnatorResultFile = 'src/test/data/cnvnator/test.cnvnator.vcf'
	static String vcgsCnvnatorFile = 'src/test/data/cnvnator/190101_A123456_0001_ML123456_TESTSAMPLE_SOMEBATCH_SOMEASSAY.cnvnator.vcf'
	

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
	public void 'Check Sample ID extracted if VCGS specific'() {
		CNVNatorResults c = new CNVNatorResults(vcgsCnvnatorFile)
		assert c.sample == 'TESTSAMPLE'
	}
	
	@Test
	public void 'Check Properties'() {
		CNVNatorResults c = new CNVNatorResults(vcgsCnvnatorFile)
		def row = c.first()
		assert row.chr == 'chrM'
		assert row.start == 1
		assert row.end == 16600
		assert row.type == 'DUP'
		assert row.sample == 'TESTSAMPLE'
		assert row.quality == 0.00809677
	}

	
	@Test
	public void "toJson Test"() {
		CNVNatorResults c = new CNVNatorResults(cnvnatorResultFile, 'someSample')
		def b = c.toListMap()
		def jsonStr = """[{"chr":"chrM","start":1,"end":16601,"sample":"someSample","quality":0.00809677,"type":"DUP"},{"chr":"chr1","start":1,"end":10001,"sample":"someSample","quality":-1.0,"type":"DEL"},{"chr":"chr1","start":10901,"end":26601,"sample":"someSample","quality":0.672392,"type":"DUP"},{"chr":"chr1","start":16201,"end":37701,"sample":"someSample","quality":0.609705,"type":"DUP"},{"chr":"chr1","start":66001,"end":133701,"sample":"someSample","quality":0.561905,"type":"DEL"},{"chr":"chr1","start":68701,"end":162101,"sample":"someSample","quality":0.833882,"type":"DEL"},{"chr":"chr1","start":130501,"end":265101,"sample":"someSample","quality":0.643947,"type":"DUP"},{"chr":"chr1","start":137301,"end":276601,"sample":"someSample","quality":0.310627,"type":"DUP"},{"chr":"chr1","start":145501,"end":298601,"sample":"someSample","quality":0.927725,"type":"DEL"},{"chr":"chr1","start":153601,"end":309701,"sample":"someSample","quality":0.986607,"type":"DEL"},{"chr":"chr1","start":177101,"end":404801,"sample":"someSample","quality":0.953488,"type":"DEL"},{"chr":"chr1","start":228001,"end":467701,"sample":"someSample","quality":0.411539,"type":"DUP"},{"chr":"chr1","start":267301,"end":585701,"sample":"someSample","quality":0.635135,"type":"DEL"},{"chr":"chr1","start":351501,"end":737501,"sample":"someSample","quality":0.954205,"type":"DEL"},{"chr":"chr1","start":426301,"end":875201,"sample":"someSample","quality":0.996778,"type":"DUP"},{"chr":"chr1","start":454201,"end":918101,"sample":"someSample","quality":0.967127,"type":"DUP"},{"chr":"chr1","start":465301,"end":936601,"sample":"someSample","quality":1.0,"type":"DUP"},{"chr":"chr1","start":471301,"end":992801,"sample":"someSample","quality":1.0,"type":"DEL"},{"chr":"chr1","start":528201,"end":1114201,"sample":"someSample","quality":0.579278,"type":"DUP"},{"chr":"chr1","start":657901,"end":1320401,"sample":"someSample","quality":0.738525,"type":"DUP"},{"chr":"chr1","start":672101,"end":1348401,"sample":"someSample","quality":0.464541,"type":"DUP"},{"chr":"chr1","start":702601,"end":1420501,"sample":"someSample","quality":0.276643,"type":"DUP"},{"chr":"chr1","start":754301,"end":1523401,"sample":"someSample","quality":0.0424754,"type":"DUP"},{"chr":"chr1","start":793401,"end":1637201,"sample":"someSample","quality":0.0206,"type":"DUP"},{"chr":"chr1","start":845101,"end":1714501,"sample":"someSample","quality":0.00220291,"type":"DUP"},{"chr":"chr1","start":870301,"end":1782701,"sample":"someSample","quality":0.0044191,"type":"DUP"},{"chr":"chr1","start":913401,"end":1925101,"sample":"someSample","quality":0.00725914,"type":"DUP"},{"chr":"chr1","start":1012801,"end":2087501,"sample":"someSample","quality":0.0222713,"type":"DUP"}]
"""
		def jsonSlurper = new JsonSlurper()
		def a = jsonSlurper.parseText(jsonStr)
		assert a.size() == b.size()
		assert b[0].sample == a[0].sample
		assert b[-1].sample == a[-1].sample
		
	}
	

}
