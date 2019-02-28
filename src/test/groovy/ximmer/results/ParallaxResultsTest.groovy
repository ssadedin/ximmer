package ximmer.results
import gngs.*

import static org.junit.Assert.*

import org.junit.Test
import groovy.json.JsonSlurper

class ParallaxResultsTest {
	
	static String pxResultFile = 'src/test/data/px/test.px.cnvs.bed'
	static String invalidPxResultFile = 'src/test/data/px/invalid_test.px.cnvs.bed'
	static String vcgsPxFile = 'src/test/data/px/190101_A123456_0001_ML123456_TESTSAMPLE_SOMEBATCH_SOMEASSAY.px.cnvs.bed'

	@Test
	public void 'TestLoadFileSuccess'() {
		ParallaxResults p = new ParallaxResults(pxResultFile).load()
	}
	
	@Test(expected = java.io.FileNotFoundException.class)
	public void 'TestLoadFileFail'() {
		ParallaxResults p = new ParallaxResults(invalidPxResultFile).load()
	}
	
	@Test
	public void 'Test ParallaxResults as a Region Object'() {
		ParallaxResults p = new ParallaxResults(pxResultFile).load()
		Region r = new Region("chr22", 16057532, 16059424)
		assert p.find { it.overlaps(r) } 
	}
	
	@Test
	public void 'Check Parallax Sample Name if not VCGS specific'() {
		ParallaxResults p = new ParallaxResults(pxResultFile).load()
		assert p.sample == 'test'
	}
	
	@Test
	public void 'Check Sample ID extracted if VCGS specific'() {
		ParallaxResults p = new ParallaxResults(vcgsPxFile).load()
		assert p.sample == 'TESTSAMPLE'	
	}
	
	@Test
	public void 'Set sample explicitly'() {
		ParallaxResults p = new ParallaxResults(pxResultFile, 'someSample').load()
		assert p.sample == 'someSample'
	}

	@Test
	public void 'Check Properties'() {
		ParallaxResults p = new ParallaxResults(pxResultFile, 'someSample').load()
		def row = p.first()
		
		assert row.type == 'DEL'
		assert row.sample == 'someSample'
		assert row.quality == 473.5896343569966
	}
	
	@Test
	public void "toJson Test"() {
		def jsonStr = """[{"chr":"chr22","start":16057532,"end":16059424,"sample":"someSample","quality":473.5896343569966,"type":"DEL"},{"chr":"chr22","start":16064346,"end":16067732,"sample":"someSample","quality":624.6264315126965,"type":"DEL"},{"chr":"chr22","start":16067971,"end":16069045,"sample":"someSample","quality":266.8606382441293,"type":"DEL"},{"chr":"chr22","start":16071061,"end":16072062,"sample":"someSample","quality":207.2024678962067,"type":"DEL"},{"chr":"chr22","start":16072251,"end":16075439,"sample":"someSample","quality":3.162423538933238E9,"type":"DEL"},{"chr":"chr22","start":16075553,"end":16077554,"sample":"someSample","quality":570.6824902981128,"type":"DEL"},{"chr":"chr22","start":16077680,"end":16079181,"sample":"someSample","quality":249.6676851419046,"type":"DEL"},{"chr":"chr22","start":16081264,"end":16083190,"sample":"someSample","quality":1.142999114981009E9,"type":"DEL"},{"chr":"chr22","start":16083570,"end":16085571,"sample":"someSample","quality":376.98063184997864,"type":"DEL"},{"chr":"chr22","start":16087854,"end":16092977,"sample":"someSample","quality":6.8532385101160194E10,"type":"DEL"},{"chr":"chr22","start":16093221,"end":16095810,"sample":"someSample","quality":6.048206382263937E10,"type":"DEL"},{"chr":"chr22","start":16101600,"end":16105262,"sample":"someSample","quality":16790.777909622462,"type":"DEL"},{"chr":"chr22","start":16107128,"end":16111370,"sample":"someSample","quality":1.4648292955042383E10,"type":"DEL"},{"chr":"chr22","start":16111454,"end":16113455,"sample":"someSample","quality":281.79784321459505,"type":"DEL"},{"chr":"chr22","start":16117487,"end":16118988,"sample":"someSample","quality":240.34669703188857,"type":"DEL"},{"chr":"chr22","start":16119456,"end":16122534,"sample":"someSample","quality":350.6047049109727,"type":"DEL"},{"chr":"chr22","start":16124371,"end":16126878,"sample":"someSample","quality":4852.639219854074,"type":"DEL"},{"chr":"chr22","start":16128822,"end":16130884,"sample":"someSample","quality":6008.02987026224,"type":"DEL"},{"chr":"chr22","start":16132260,"end":16133261,"sample":"someSample","quality":577.7840329093116,"type":"DEL"},{"chr":"chr22","start":16133381,"end":16135279,"sample":"someSample","quality":270.5669440528906,"type":"DEL"}]
"""
		def jsonSlurper = new JsonSlurper()
		def a = jsonSlurper.parseText(jsonStr)
		ParallaxResults p = new ParallaxResults(pxResultFile, 'someSample').load()
		def b = p.toListMap()
		assert a.size() == b.size()
		assert b[0].sample == a[0].sample
		assert b[-1].sample == a[-1].sample
	}
}
