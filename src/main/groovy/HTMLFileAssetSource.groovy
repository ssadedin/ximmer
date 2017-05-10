class HTMLFileAssetSource extends HTMLAssetSource {
    
    File templateParentDir
    
    HTMLFileAssetSource(File templateParentDir) {
        this.templateParentDir = templateParentDir
    }
    
    @Override
    public String getText(String fileName) {
      new File(templateParentDir,fileName).newInputStream().withStream { it.text }
    }
}
