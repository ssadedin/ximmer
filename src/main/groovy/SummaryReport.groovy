import groovy.text.SimpleTemplateEngine
import groovy.util.logging.Log

import java.nio.file.Files
import java.nio.file.StandardCopyOption

import gngs.*

@Log
class SummaryReport {
    
    Ximmer ximmer
    
    File dir
    
    /**
     * These assets are copied to the same folder as the HTML report
     */
    static List<String> SUMMARY_HTML_ASSETS = [
        'jquery-2.1.0.min.js',
        'summary_report.js',
        'DOMBuilder.dom.min.js',
        'jquery-ui.min.js',
        'cnv_diagram.js',
        'require.js',
        'cnv_report.js',
        'jquery-ui.css',
        'nv.d3.css',
        'c3.css',
        'N3Components.min.css',
        'ximmer.css',
        'cnv_report.css'
    ]
    
    /**
     * These assets are loaded via requirejs and thus should NOT
     * have script tags written
     */
    static List<String> NO_TAG_HTML_ASSETS = [
        'summary_report.js',
        'N3Components.min.js',
        'roccurve.js',
        'require.js',
    ]
    
    SummaryReport(Ximmer ximmer) {
        this.ximmer = ximmer
        this.dir = ximmer.outputDirectory
    }
    
    /**
     * Writes the top level summary report containing simulation information,
     * aggregate statistics, as well as tabs to view the lower level CNV reports.
     * 
     * @param analysis
     */
    void write(AnalysisConfig analysis) {
        
        String summaryHTML = generateSummary(analysis)
        
        String analysisName = analysis.analysisName
        
        List runDirectories = ximmer.runs*.value*.runDirectory;
        
        writeEncodedCallerReports(runDirectories, analysisName)
        
        log.info("Generating HTML Report ...")
        File mainTemplate = new File("$ximmer.ximmerBase/src/main/resources/index.html")
        
        String outputName = analysisName + ".html"
        
        HTMLAssetSource source = new HTMLClassloaderAssetSource()
        HTMLAssets assets = new HTMLAssets(source, dir)
        for(asset in SUMMARY_HTML_ASSETS.grep { !(it in NO_TAG_HTML_ASSETS) })
            assets << new HTMLAsset(source:asset)
            
        String assetPayload = assets.render()
        
        copyResources()
        
        new File(dir, outputName).withWriter { w ->
            SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
            templateEngine.createTemplate(mainTemplate.newReader()).make(
                analysisName : analysisName,
                runDirectories: runDirectories,
                outputDirectory : dir.name,
                summaryHTML : summaryHTML,
                analysisConfig: analysis,
                callers: ximmer.callerIds,
                callerIdMap: ximmer.callerIdMap,
                enableTruePositives: ximmer.enableTruePositives,
                assets: assetPayload,
                simulation_type: ximmer.cfg.simulation_type,
                config: ximmer.cfg
            ).writeTo(w)
        }
    }
    
    
    /**
     * Write out each CNV report, encoded as a javascript string using base 64.
     * <p>
     * This allows the main ximmer report to load and inject these reports as HTML on the fly
     * without having all the HTML for each report (which can be huge) loaded into memory.
     * 
     * @param runDirectories
     * @param analysisName
     */
    void writeEncodedCallerReports(List<File> runDirectories, String analysisName) {
        runDirectories.eachWithIndex { File runDir, runIndex  -> 
            File b64File = new File(dir, runDir.name+'/' + analysisName +'/report/cnv_report.b64.js')
            log.info "Writing base 64 encoded CNV report to " + b64File
            b64File.withWriter { w ->
                w.println 'var cnvReportHTML = "' + new File(dir, runDir.name+'/' + analysisName +'/report/cnv_report.html').text.bytes.encodeBase64() + '";'
            }
        }
    }
    
    /**
     * Returns HTML for the summary / overview tab
     * 
     * @return
     */
    String generateSummary(AnalysisConfig analysisCfg) {
        
        List<Region> cnvs = ximmer.readCnvs()
        
        // Load the results
        List<RangedData> results = ximmer.runs*.value*.runDirectory.collect { runDir ->
            new RangedData(new File(runDir,"$analysisCfg.analysisName/report/cnv_report.tsv").path).load([:], { r ->
                    (ximmer.callerIds + ["truth"]).each { r[it] = (r[it] == "TRUE") }
          })
        }
        
        plotCNVSizeHistograms(cnvs)
        
//        List<IntRange> sizeBins = [0..<200, 200..<500, 500..<1000, 1000..<2000,2000..<10000]
//        Map<IntRange,Integer> binnedCounts = 
        
        Map<String,Integer> callerCounts = ximmer.callerIds.collectEntries { caller -> 
            [caller, results.sum { r -> r.count { cnv -> cnv.truth && cnv[caller] } }] 
        }
        
        File mainTemplate = new File("$ximmer.ximmerBase/src/main/resources/summary.html")
        
        StringWriter result = new StringWriter()
        SimpleTemplateEngine templateEngine = new SimpleTemplateEngine()
        templateEngine.createTemplate(mainTemplate.newReader()).make(
            cnvs: cnvs,
            bamFiles: ximmer.bamFiles,
            batch_name : dir.name,
            callers: ximmer.callerIds,
            simulation_type: ximmer.cfg.simulation_type,
            results: results,
            callerCounts: callerCounts
        ).writeTo(result)
        return result.toString()
    }
    
    void plotCNVSizeHistograms(List<Region> cnvs) {
        
        if(!ximmer.enableTruePositives)
            return
        
        File combinedCnvs = writeCombinedCNVInfo(cnvs)
        
        ximmer.runPython([:], dir, new File("$ximmer.ximmerBase/src/main/python/cnv_size_histogram.py"), 
                  [combinedCnvs.absolutePath, new File(dir,"cnv_size_histogram.png").absolutePath])
    }
    
    File writeCombinedCNVInfo(List<Region> cnvs) {
        
       if(!dir.exists())
           dir.mkdirs()
       
       File cnvFile = new File(dir,"combined_cnvs.tsv")
       cnvFile.withWriter { w ->
            w.println([ "chr","start","end","sample","tgbp","targets" ].join("\t"))
            
            if(!ximmer.enableTruePositives) 
                return   
            
            cnvs.collect { 
                [
                    it.chr, 
                    it.from, 
                    it.to,
                    it.sample, 
                    it.targeted,
                    it.targets
                ]   
            }.each {
                w.println it.join("\t")
            }
        }
       return cnvFile
    }
    
    
    
    void copyResources() {
        for(String asset in SUMMARY_HTML_ASSETS) {
            File assetFile = new File(dir,asset)
            if(!assetFile.exists()) {
                File sourceFile = new File("$ximmer.ximmerBase/src/main/resources/$asset")
                log.info "Copy $sourceFile => $assetFile"
                Files.copy(sourceFile.toPath(), 
                           assetFile.toPath())
            }
        }
        
        // Automatically copy everything in src/lib/js 
        List<File> jsFiles = new File("$ximmer.ximmerBase/src/main/js").listFiles().grep { it.name.endsWith('.js') }
        for(File sourceFile in jsFiles) {
            File assetFile = new File(dir,sourceFile.name)
            log.info "Copy $sourceFile => $assetFile"
            Files.copy(sourceFile.toPath(), assetFile.toPath(), StandardCopyOption.REPLACE_EXISTING)
        }
    }
    

}
