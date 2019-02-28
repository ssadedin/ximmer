package ximmer.results

import static org.junit.Assert.*

import groovy.json.JsonSlurper
import org.junit.Test

class CodexResultsTest {
	
	static String cdxResultFile = 'src/test/data/cdx/test.cdx.tsv'

	@Test
	public void 'TestLoadFileSuccess'() {
		CodexResults c = new CodexResults(cdxResultFile).load()
		
	}
	
	@Test(expected = java.io.FileNotFoundException.class)
	public void 'TestLoadFileFail'() {
		CodexResults c = new CodexResults(cdxResultFile + 'blah').load()
	}
	
	@Test
	public void 'Check Properties'() {
		CodexResults c = new CodexResults(cdxResultFile).load()
		def row = c.first()
		assert row.sample == '10X000626-XZ'
		assert row.chr == 'chr1'
		assert row.type == 'DEL'
		assert row.quality == 14.564
	}
	
	@Test
	public void "toJson Test"() {
		String jsonStr = """[{"chr":"chr1","start":78383809,"end":78390914,"sample":"10X000626-XZ","quality":14.564,"type":"DEL"},{"chr":"chr1","start":112319647,"end":112321114,"sample":"10X000507-OC","quality":23.982,"type":"DEL"},{"chr":"chr1","start":156084460,"end":156100564,"sample":"10X001329SB","quality":166.756,"type":"DEL"},{"chr":"chr10","start":69966526,"end":69971773,"sample":"10X000826TG","quality":169.764,"type":"DUP"},{"chr":"chr10","start":92677488,"end":92679025,"sample":"10X000025JK","quality":12.393,"type":"DEL"},{"chr":"chr10","start":92678621,"end":92679025,"sample":"10X001329SB","quality":4.657,"type":"DEL"},{"chr":"chr14","start":90870722,"end":90874619,"sample":"10X000507-OC","quality":12.836,"type":"DEL"},{"chr":"chr17","start":12620662,"end":12623778,"sample":"10X001402CM","quality":18.224,"type":"DEL"},{"chr":"chr17","start":28295808,"end":28296385,"sample":"10X000647-DG","quality":8.095,"type":"DEL"},{"chr":"chr18","start":3173935,"end":3176132,"sample":"10X000499-SC","quality":16.671,"type":"DEL"},{"chr":"chr19","start":16594740,"end":16596068,"sample":"10X001658BoEP","quality":16.346,"type":"DEL"},{"chr":"chr2","start":179549632,"end":179550057,"sample":"10X001329SB","quality":21.818,"type":"DEL"},{"chr":"chr2","start":179549056,"end":179549476,"sample":"10X001329SB","quality":19.637,"type":"DUP"},{"chr":"chr2","start":179539040,"end":179539840,"sample":"10X001329SB","quality":18.237,"type":"DEL"},{"chr":"chr2","start":189853315,"end":189854175,"sample":"10X000499-SC","quality":12.754,"type":"DUP"},{"chr":"chr4","start":113825610,"end":113970968,"sample":"10X000507-OC","quality":10.843,"type":"DEL"},{"chr":"chr6","start":123696749,"end":123759267,"sample":"10X000046EB","quality":10.416,"type":"DEL"},{"chr":"chr7","start":81596911,"end":81598289,"sample":"10X000626-XZ","quality":11.713,"type":"DUP"},{"chr":"chrX","start":70360445,"end":70361269,"sample":"10X000626-XZ","quality":45.932,"type":"DEL"},{"chr":"chrX","start":70360445,"end":70361269,"sample":"10X000507-OC","quality":50.93,"type":"DEL"},{"chr":"chrX","start":70360445,"end":70361269,"sample":"10X000647-DG","quality":64.303,"type":"DEL"}]
"""
		def jsonSlurper = new JsonSlurper()
		def a =  jsonSlurper.parseText(jsonStr)
		
		CodexResults c = new CodexResults(cdxResultFile).load()
		def codexJsonStr = c.toJson()
		def b = c.toListMap()
		assert a.size() == b.size()
		assert b[0].sample == a[0].sample
		assert b[-1].sample == a[-1].sample
	}
}
