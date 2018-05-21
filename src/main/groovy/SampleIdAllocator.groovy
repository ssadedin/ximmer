import groovy.util.logging.Log

@Singleton
@Log
class SampleIdAllocator {
    
    int offset = 1
    
    String prefix = 'XS'
    
    Map<String,String> allocatedIds = [:]
    
    boolean anonymise = false
    
    synchronized String newSampleId(String sampleId) {
        
        if(allocatedIds[sampleId])
            return allocatedIds[sampleId]
        
        String newId = anonymise ? String.format(prefix + '%03d', offset++) : sampleId
        allocatedIds[sampleId] = newId
        return newId
    }

}
