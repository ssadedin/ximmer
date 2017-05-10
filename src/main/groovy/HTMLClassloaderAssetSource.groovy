
class HTMLClassloaderAssetSource extends HTMLAssetSource {
    @Override
    public String getText(String fileName) {
        InputStream stream = getClass().classLoader.getResourceAsStream(fileName)
        
        if(stream == null)
            throw new FileNotFoundException("Cannot locate file $fileName in classpath")
            
        stream.withStream { it.text }
    }
}



