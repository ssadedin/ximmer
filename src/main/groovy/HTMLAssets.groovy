import java.io.InputStream
import java.util.logging.Logger

import groovy.transform.CompileStatic
import groovy.util.logging.Log


/**
 * Outputs a script tag and also writes a file it refers to into the output
 * 
 * @author Simon Sadedin
 */
@Log
class HTMLAssets {
    
    HTMLAssetSource source;
    
    List<HTMLAsset> assets = []
    
    File outputDir
    
    boolean inline = false
    
    HTMLAssets(HTMLAssetSource source, File outputDir) {
        this.source = source;
        this.outputDir = outputDir
    }

    @CompileStatic
    HTMLAssets add(HTMLAsset asset) {
        if(asset.name == null)
            asset.name = asset.source
        assets << asset
        return this
    }
    
    HTMLAssets leftShift(HTMLAsset asset) {
        add(asset)
    }
    
    @CompileStatic
    String render() {
        Logger lg = log
        assets.collect { HTMLAsset script ->
            String code = source.getText(script.source)
            if(inline) {
                if(script.name.endsWith('.js')) {
                    """<script type="text/javascript">\n${code}\n</script>\n"""
                }
                else {
                    """<style type="text/css">\n${code}\n</style>\n"""
                }
            }
            else {
                File assetFile = new File(outputDir, script.name)
                lg.info "Copy $script => $assetFile"
                assetFile.text = code
                if(script.name.endsWith('.js')) {
                    """<script type="text/javascript" src='$script.name'></script>\n"""
                }
                else { // assume stylesheet
                    """<link rel="stylesheet" href='$script.name'>"""
                }
            }
        }.join("\n")
    }
}
