import java.util.regex.Pattern

import com.github.scribejava.core.builder.api.DefaultApi10a

import gngs.*
import groovy.json.JsonSlurper
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord

/**
 * Utility to send raw JSON from stdin to CXP URL, with appropriate signing
 * <p>
 * Expects to find OAuth credentials in a file called <code>.cxp/credentials</code>
 * which should have the form:
 * 
 * <pre>
 * apiKey:  anapikey
 * apiSecret: theapisecret
 * accessToken: anaccesstoken
 * accessSecret: theaccesssecret
 * </pre>
 * @author Simon Sadedin
 */
@Log
class RawPostToCXP extends ToolBase {

    private WebService ws
    
    @Override
    public void run() {

        log.info "Posting to $opts.url"
           
        ws = new WebService(opts.url)
        ws.autoSlash = false
        ws.credentialsPath = ".cxp/credentials"
        
        def data = new JsonSlurper().parse(System.in)
        
        log.info 'Posting data: ' + data
        
        ws.post(data)
    }
    
    static void main(String [] args) {
        cli('RawPostToCXP -url <full post url>', args) {
            url 'URL to post to', args:1, required: true
        }
    }
}
