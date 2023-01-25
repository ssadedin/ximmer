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
class PostToCXPWGS extends ToolBase {
    private WebService ws
    
    private Ximmer ximmer
    
    private Map<String,String> sampleSexes = [:]
    
    private Map<String,List<String>> bamFiles = [:]

    private Map<String,List<String>> vcfFiles = [:]
    
    private ConfigObject cfg
    
    private WebService createBamService

    private WebService createVcfService
    
    static void main(String [] args) {
        cli('PostToCXPWGS -cxp <CXP URL> -analysis <analysis directory> -sex <sample_id:SEX> -bam <sample_id:bam_file>', args) {
            project 'The project GUID to associate the analyses to', args: 1, required: true
            analysis 'The zip file of the analysis to import', args:1, required: true
            sex 'Provide sample_id:SEX', args:UNLIMITED, required: true
            bam 'Provide sample_id:BAM', args:UNLIMITED, required: true
            qc 'Zip file containing QC files to import', args:1, required: true
            cxp 'CXP Url', args:1, required: true
            target 'Provide the path to target region', args:1, required: true
            batch 'Fullpath to batch', args:1, required: false
            test 'Test mode: do not actually send requests, just print them out'
            vcf 'Provide sample_id:VCF', args: UNLIMITED, required: false
        }
    }

    public void getSampleFilesFromColonArgs(args, typeName, sampleFileMap, Closure filePredicate = null) {
        if (!args) return
        // Sample Files common in from args have the particular structure sampleId:[Fullpath]
        for(arg in args) {
            if (!arg.contains(':')) {
                throw new IllegalArgumentException("Please provide -${typeName} sample_id:${typeName.toUpperCase()}")
            } 
            def parts = arg.tokenize(':')
            def sampleId = parts[0]
            def fullpath = parts[1]
            if (!new File(fullpath).exists()) {
                throw new IllegalArgumentException("${typeName}file: [$fullpath] doesn't exist")
            }
            if (filePredicate)
                filePredicate(fullpath)
            // Only accept if sampleId has already been inserted into bamFiles
            if (!sampleFileMap.keySet().contains(sampleId)) 
                throw new IllegalArgumentException("${typeName}file with sample: $sampleId provided without sex")
            log.info "Adding to ${typeName}files: [sample:$sampleId] [${typeName}File:$fullpath]"
            sampleFileMap[sampleId] << new File(fullpath).absolutePath
        }
    }

    @Override
    public void run() {
        log.info "Target Region: ${opts.target}"
        String assay = new File(opts.target).name.replaceAll('.bed$', '')
        log.info "Assay: $assay"
        log.info "Analysis zip: ${opts.analysis}"
        log.info "Target project: ${opts.project}"

        String projectGuid = opts.project

        File batchDir;
        if (opts.batch) {
            batchDir = new File(opts.batch).absoluteFile
        } else {
            batchDir = new File('.').absoluteFile.parentFile
        }
        log.info "batchDir: ${batchDir}"
        
        if(!opts.cxp.contains('http')) throw new IllegalArgumentException("Invalid CXP URL used")
        log.info "CXP url: ${opts.cxp}"
        
        if (opts.sexs) {
            for(sexArg in opts.sexs) {
                if (!sexArg.contains(':')) {
                    throw new IllegalArgumentException("Please provide -sex sample_id:SEX")
                } 
                def parts = sexArg.tokenize(':')
                def sampleId = parts[0]
                def sex = parts[1].toUpperCase()
                if (!['MALE','FEMALE','UNKNOWN'].contains(sex)) 
                    throw new IllegalArgumentException("Sex provided was not MALE|FEMALE|UNKNOWN")
                log.info "Adding to sampleSexes: [sample:$sampleId] [sex:$sex]"
                this.sampleSexes[sampleId] = sex
                this.bamFiles[sampleId] = []
                this.vcfFiles[sampleId] = []
            }

            getSampleFilesFromColonArgs(opts.bams, 'bam', this.bamFiles) { bamFile -> 
                if (!new File(bamFile + ".bai").exists()) {
                   throw new IllegalArgumentException("Bam Index file for: [$bamFile] doesn't exist")
               }
            }

            getSampleFilesFromColonArgs(opts.vcfs, 'vcf', this.vcfFiles)

            this.bamFiles.each { k, v -> if(!v) throw new IllegalArgumentException("Couldn't find any bamFiles for sample: $k")}
        } 
        
        ws = new WebService(opts.cxp)
        ws.autoSlash = true
        ws.credentialsPath = ".cxp/credentials"
        this.createBamService = (ws / 'dataasset/create/bam/')
        this.createVcfService = (ws / 'dataasset/create/vcf/')
        
        this.registerBAMFiles(batchDir, assay)
        this.registerVCFs(batchDir, assay)
        this.postAnalysis(batchDir, assay, projectGuid)
        
    }
    
    /**
     * Post details of the analysis (CNV calls and QC results) to the CXP API end point
     * 
     * @param batchDir
     * @param assay
     */
    void postAnalysis(File batchDir, String assay, String projectGuid) {
        String sequencer = new SAM(new File(this.bamFiles.values()[0][0])).withIterator { i -> 
            SAMRecord r = i.next()
            return r.readName.tokenize(':')[0].stripMargin('@')
        }
        
        log.info "Found sequencer: $sequencer"
        
        String batchIdentifier = batchDir.name
        // An analysis needs a batch, so create one?
        List batch = (ws / 'batch').get(identifier:batchIdentifier)
        if(batch) {
            println "Found batch $batch"
        }
        else {
            log.info "Creating new batch $batchIdentifier"
            def data =  [
               metadata: [:],
               identifier: batchIdentifier
            ]

            if(opts.test) {
                log.info "Would post:\n\n$data\n\n to ${ws/'batch'}"
            }
            else {
                batch = [(ws / 'batch').post(data)]
            }
        }
        
       
        Map data = [
            identifier: batchDir.absolutePath,
            project_guid: projectGuid,
            assay: assay,
            sequencer: sequencer,
            // samples: this.bamFiles.keySet(),
            samples: this.bamFiles.collectEntries { [it.key, it.value[0] ] },
            batch_id: batch[0].id,
            results: new File(opts.analysis).absolutePath,
            control_samples: [],
            analysis_samples: this.bamFiles.keySet(),
            qc: new File(opts.qc).absolutePath
        ]
        
        log.info "Sending to analysis/import data: $data"
       
        WebService importService = ws / 'analysis/import'
        
        if(opts.test) {
            log.info "Would post\n\n$data to \n\n$importService.endPoint"
        } 
        else {
            importService.post(data) 
        }
        
    }
    
    /**
     * For each BAM file in the Ximmer analysis, register it with CXP
     * 
     * @param batchDir
     * @param assay
     */
    void registerBAMFiles(File batchDir, String assay) {
        bamFiles.each { sampleId, bams -> 
            bams.each { 
                bamFile -> registerBAM(sampleId, bamFile, assay, batchDir) 
            }
        }
    }


    /**
     * For each VCF file in the Ximmer analysis, register it with CXP
     * 
     * @param batchDir
     * @param assay
     */
    void registerVCFs(File batchDir, String assay) {
        vcfFiles.each { sampleId, vcfs -> 
            vcfs.each { 
                vcfFile -> registerVCF(sampleId, vcfFile, assay, batchDir) 
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
    void registerBAM(String sampleId, String bamFile, String assay, File batchDir) {
        //log.info "registerBAM with: sampleId:$sampleId, bamFile:$bamFile, assay:$assay, batchDir:$batchDir"
      
        Map data = [
            'fullpath': new File(bamFile).absolutePath,
            'sample': sampleId,
            'filetype': '.bam',
            'alt_sample_id': sampleId,
            'sex': sampleSexes[sampleId],
            'batch': batchDir.name,
            'sequencer': 'Test',
            'assay': assay,
            'batch_date': batchDate(batchDir.lastModified()),
            'metadata': [:]
        ]
        
        println data
        
        
        if(opts.test) {
            log.info "Would post $data to $createBamService.endPoint"
        }
        else {
            createBamService.post(data)
        }
    }

    /**
     * Register the given VCF file with CXP
     * <p>
     * NOTE: sex is inferred from the VCF file itself
     *
     * @param sampleId
     * @param vcfFile
     * @param assay
     * @param batchDir
     */
     void registerVCF(String sampleId, String vcfFile, String assay, File batchDir) {
        Map data = [
            'fullpath': new File(vcfFile).absolutePath,
            'sample': sampleId,
            'filetype': 'vcf',
            'sex': sampleSexes[sampleId],
            'batch': batchDir.name,
            'assay': assay,
        ]
        
        println data
        
        if(opts.test) {
            log.info "Would post $data to $createVcfService.endPoint"
        }
        else {
            createVcfService.post(data)
        }
    }

    
    
    String batchDate(long timeMs) {
       new Date(timeMs).format('YYYY-MM-dd') 
    }


}
