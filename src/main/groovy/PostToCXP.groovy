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
class PostToCXP extends ToolBase {
    private WebService ws
    
    private Ximmer ximmer
    
    private Map<String,String> sampleSexes = [:]
    
    private ConfigObject cfg
    
    private WebService createBamService

    @Override
    public void run() {
        log.info "Loading configuration from ${new File(opts.c).absolutePath}"
        cfg = new ConfigSlurper().parse(new File(opts.c).text)
        log.info "Configuration parsed."
       
        ximmer = new Ximmer(cfg, "cnv", false)
        ximmer.initialiseRuns()
            
        ws = new WebService(opts.cxp)
        ws.autoSlash = true
        ws.credentialsPath = ".cxp/credentials"
        
        this.createBamService = (ws / 'dataasset/create/bam/')
        
        File batchDir = new File('.').absoluteFile.parentFile
        String assay = opts.assay?:new File(cfg.target_regions).name.replaceAll('.bed$','')
        
        // Step 1: Infer the sexes
        this.inferSexes(ximmer.bamFiles*.value)
        
        log.info "Inferred sexes: \n" + Utils.table(["Sample","Sex"], [ sampleSexes*.key,sampleSexes*.value].transpose())
        
        // Step 2: Make sure the BAM files are registered
        this.registerBAMFiles(assay)
        
        // Step 3: Register a new analysis - or can this happen automatically for a result that doesn't have an analysis?
        this.postAnalysis(batchDir, assay)
    }
    
    /**
     * Post details of the analysis (CNV calls and QC results) to the CXP API end point
     * 
     * @param batchDir
     * @param assay
     */
    void postAnalysis(File batchDir, String assay) {
        
        String sequencer = ximmer.bamFiles*.value[0].withIterator { i -> 
            SAMRecord r = i.next()
            return r.readName.tokenize(':')[0].stripMargin('@')
        }
        
        String batchIdentifier = opts.batch?:batchDir.name
        
        // An analysis needs a batch, so create one?
        List batch = (ws / 'batch').get(identifier:batchIdentifier)
        if(batch) {
            println "Found batch $batch"
        }
        else {
            log.info "Creating new batch $batchIdentifier"
            batch = [(ws / 'batch').post(
                metadata: [:],
                identifier: batchIdentifier,
                date: batchDate(batchDir.lastModified())
            )]
        }
        
        List samplesToSubmit = ximmer.bamFiles*.key
        def requiredSex = this.cfg.getOrDefault('filter_to_sex',false)
        if(requiredSex) {
            samplesToSubmit = samplesToSubmit.grep { sampleSexes[it] == requiredSex }
        }
        
        Map data = [
            identifier: batchDir.absolutePath,
            assay: assay,
            sequencer: sequencer,
//            samples: ximmer.bamFiles, // samplesToSubmit,
            samples: samplesToSubmit.collectEntries { [ it, ximmer.bamFiles[it] ] },
            batch_id: batch[0].id,
            results: new File(opts.analysis).absolutePath,
            control_samples: [],
            analysis_samples: ximmer.bamFiles*.key,
            qc: new File(opts.qc).absolutePath
        ]
        
        WebService importService = ws / 'analysis/import'
        
        if(opts.test) {
            log.info "Would post $data to $createBamService.endPoint"
        } 
        else {
            importService.post(data) 
        }
        
    }
    
    /**
     * For each BAM file in the Ximmer analysis, register it with CXP
     * 
     * @param assay
     */
    void registerBAMFiles(String assay) {
        for(SAM bam in ximmer.bamFiles*.value) {
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
        Regions karyoRegions = ximmer.targetRegion.grep { it.chr.replaceAll(chrStart,'') in karyoChr } as Regions
        karyoRegions = karyoRegions.thin(200, 50)
        if(!opts.sex)
            log.info "Karyotyping using ${Utils.humanBp(karyoRegions.size())} consisting of ${karyoRegions.numberOfRanges} target regions from total ${Utils.humanBp(ximmer.targetRegion.size())} in target region"
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
        cli('PostToCXP -c <Ximmer Config> -cxp <CXP URL> <analysis directory>', args) {
            c 'Ximmer Configuration File', args:1, required: true
            analysis 'The zip file of the analysis to import', args:1, required: true
            qc 'The zip file containing QC files to import', args:1, required: true
            cxp 'Base URL to CXP server', args:1, required: true
            assay 'The assay to register the samples under (default: name of BED file)', args:1, required: false
            batch 'CNV calling batch identifier (default: name of current directory', args:1, required: false
            test 'Do not actually post data, just show what would be posted', required: false
            sex 'Specify sex for all samples', args: 1, required: false
            bambatch 'The batch name associated to the BAM files', args:1, required: false
            sequencer 'The sequencer used to create the raw data in the BAM files', args:1, required: false
        }
    }
}
