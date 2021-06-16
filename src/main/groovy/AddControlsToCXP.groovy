import java.util.regex.Pattern

import com.github.scribejava.core.builder.api.DefaultApi10a

import gngs.*
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord

/**
 * Utility to load analysis results from Ximmer into a CXP API end point
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
class AddControlsToCXP extends ToolBase {

    private WebService ws
    
    private Map<String,String> sampleSexes = [:]
    
    private WebService createBamService
    
    private List<SAM> bamFiles
    
    private Regions targetRegions

    @Override
    public void run() {
        
        ws = new WebService(opts.cxp)
        
        // necessary to ensure URLs are normalised such that they pass signature verification (Oauth)
        ws.autoSlash = true
        ws.credentialsPath = ".cxp/credentials"
        
        List<Map> assays = ws.path('assay').get(identifier:opts.assay)
        if(!assays) 
            throw new IllegalArgumentException("The assay $opts.assay could not be resolved by CXP server at $opts.cxp")
            
        this.targetRegions = new BED(assays[0].fullpath).load()

        this.createBamService = (ws / 'dataasset/create/bam/')
        
        // Get the assay target region bed file
        File batchDir = new File('.').absoluteFile.parentFile
        
        this.bamFiles = opts.arguments().collect { new SAM(it) }

        // Step 1: Infer the sexes
        this.inferSexes(this.bamFiles)

        log.info "Inferred sexes: \n" + Utils.table(["Sample","Sex"], [ sampleSexes*.key,sampleSexes*.value].transpose())
        
        // Step 2: Make sure the BAM files are registered
        this.registerBAMFiles(opts.assay)
    }
    
    /**
     * For each BAM file in the Ximmer analysis, register it with CXP
     * 
     * @param assay
     */
    void registerBAMFiles(String assay) {
        for(SAM bam in this.bamFiles) {
           registerBAM(bam, assay)
        }
    }
    
    /**
     * Infer sex for each BAM file. Primarily just so we don't have to worry about them being
     * declared.
     */
    void inferSexes(List<SAM> bamFiles) {

        log.info "Inferring sexes for ${bamFiles.size()} bam files"
        List karyoChr = ['1','X','Y']
        Pattern chrStart = ~'^chr'
        Regions karyoRegions = this.targetRegions.grep { it.chr.replaceAll(chrStart,'') in karyoChr } as Regions
        if(karyoRegions.numberOfRanges > 300)
            karyoRegions = karyoRegions.thin(200, 50)
        if(!opts.sex)
            log.info "Karyotyping using ${Utils.humanBp(karyoRegions.size())} consisting of ${karyoRegions.numberOfRanges} target regions from total ${Utils.humanBp(targetRegions.size())} in target region"
        else
            log.info "Sex pre-specified as $opts.sex for all samples"
        
        for(bam in bamFiles) {
            String sampleId = bam.samples[0]
            if(opts.sex) {
                sampleSexes[sampleId] = opts.sex
            }
            else {
                SexKaryotyper karyotyper = new SexKaryotyper(bam, karyoRegions)
                karyotyper.run()
                log.info "Sample $sampleId => $karyotyper.sex"
                sampleSexes[sampleId] = karyotyper.sex.toString()        
            }
        }
    }
    
    /**
     * Register the given BAM file with CXP
     * <p>
     * NOTE: sex is inferred from the BAM file itself. 
     *
     * @param bam
     * @param assay
     */
    void registerBAM(SAM bam, String assay) {
        
        File bamDir = bam.samFile.absoluteFile.parentFile
        File bamBatchDir = bamDir.parentFile
        String sampleId = bam.samples[0]
        
        Map data = [
            'fullpath': bam.samFile.absolutePath,
            'sample': sampleId,
            'filetype': '.bam',
            'alt_sample_id': sampleId,
            'sex': sampleSexes[sampleId],
            'batch': opts.bambatch?:(bamBatchDir.name + '_' + bamDir.name),
            'sequencer': opts.sequencer?:'Unknown',
            'assay': assay,
            'batch_date': batchDate(bamBatchDir.lastModified()),
            'metadata': [:]
        ]
        
        if(opts.test) {
            log.info "Would post $data to $createBamService.endPoint"
        }
        else {
            createBamService.post(data)
        }
    }
    
    String batchDate(long timeMs) {
       new Date(timeMs).format('YYYY-MM-dd') 
    }

    static void main(String [] args) {
        cli('AddControlsToCXP -cxp <CXP URL> <bam files>', args) {
            cxp 'Base URL to CXP server', args:1, required: true
            assay 'The assay to register the samples under (must be registered already in CXP)', args:1, required: true
            test 'Do not actually post data, just show what would be posted', required: false
            sex 'Specify sex for all samples', args: 1, required: false
            batch 'The batch name associated to the BAM files', args:1, required: false
            sequencer 'The sequencer used to create the raw data in the BAM files', args:1, required: false
        }
    }
}
