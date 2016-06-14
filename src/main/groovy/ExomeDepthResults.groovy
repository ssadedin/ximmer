
class ExomeDepthResults extends CNVResults {

    public ExomeDepthResults(String sourceFile) {
        super(sourceFile, 6, 4, 5);
        this.separator = "\t"
    }
    
    ExomeDepthResults load(Map options=[:], Closure c=null) {
        super.load(options+[separator:"\t"]) { r ->
            r.quality = r.BF
            r.type = (r.type == "deletion") ? "DEL" : "DUP"
        }
    }
    
    public static void main(String [] args) {
        def ed = new ExomeDepthResults("testdata/asdsibs.exome_depth.cnvs.tsv").load()
        println ed[0].sample + " with quality " + ed[0].quality
    }
}
