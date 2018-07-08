package ximmer

import gngs.ToolBase
import gngs.XPos
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.sql.Sql
import groovy.util.logging.Log

@Log
class CNVDB extends ToolBase {
    
    static Map<Integer, Map<String,List<String>>>  schema = [
        1 : [
            """
                create table schema_meta_data (
                       id integer primary key asc, 
                       schema_version integer NOT NULL
                );    
            """,
            """
              insert into schema_meta_data values(1,1);
            """,
            """
                create table cnv ( 
                        id integer primary key asc, 
                        sample_id text NOT NULL,
                        chr text NOT NULL,
                        start integer,            -- the start of actual DNA change
                        end integer,              -- end end of actual DNA change
                        startx integer NOT NULL,  -- xpos of start
                        endx integer NOT NULL,     -- xpos of end
                        metadata text
                );
            """,
            """
                CREATE INDEX start_idx ON cnv(startx);
            """,
            """
                CREATE INDEX end_idx ON cnv(endx);
            """,
            """
                CREATE INDEX sample_idx ON cnv (sample_id);
            """ 
            ]
    ]
    
    
    Sql db

    @Override
    public void run() {
        
        log.info "Creating / checking database $opts.db"
        db = createDb(opts.db)
        
        log.info "Reading cnvs from ${opts.arguments().join(', ')}"
        List<Map> cnvs = opts.arguments().collect { 
            new JsonSlurper().parseText(new File(it).text)
        }.flatten()
        
        log.info "Loaded ${cnvs.size()} CNV calls"
        cnvs.each {
            insertCNV(it)
        }
        
        db.commit()
    }
    
    void insertCNV(Map cnv) {
        db.execute """
            insert into cnv(id, sample_id, chr, start, end, startx, endx, metadata)
                    values (NULL, $cnv.sample, $cnv.chr, $cnv.start, $cnv.end, ${XPos.computePos(cnv.chr, cnv.start)}, ${XPos.computePos(cnv.chr, cnv.end)}, ${JsonOutput.toJson(cnv)})
        """
    }
    
    Sql createDb(String fileName) {
        Class.forName("org.sqlite.JDBC")
        Sql db = Sql.newInstance("jdbc:sqlite:$fileName")
        
        db.cacheConnection = true
        
        // For performance
        db.execute("PRAGMA synchronous = 0;")
        db.execute("PRAGMA journal_mode = WAL;") 
        db.connection.autoCommit = false
        
        checkSchema(db)
        
        return db
    }
    
    /**
     * Execute the given statements to upgrade the schema
     * 
     * @param statements
     */
    void upgradeSchema(Sql db, int toVersion) {
        log.info "Schema is out of date: upgrading to version $toVersion"
        db.execute("BEGIN TRANSACTION;")
        List<String> statements = schema[toVersion] 
        for(String s in statements) {
            log.info "Executing upgrade statement: $s"
            db.execute(s)
        }
        db.execute("update schema_meta_data set schema_version=$toVersion")
        log.info "Committing changes"
        db.execute("COMMIT;")
    } 
    
    def schemaInfo = null
    
    /**
     * Iterate over the schema versions and check if they exist
     */
    void checkSchema(Sql db) {
        
        // Get the current schema information
        try {
            schemaInfo = db.firstRow("select * from schema_meta_data")
        }
        catch(Exception e) {
            schemaInfo = [id:0, schema_version:0] // 0 meaning, "does not exist"
        }
        
        List schemaVersions = (schema.keySet() as List).sort()
        schemaVersions.each { version ->
            if(schemaInfo.schema_version < version)  {
                upgradeSchema(db, version)
                schemaInfo.schema_version = version
            }
        } 
    }
    
    static void main(String [] args) {
        cli('CNVDB <CNV json / js', args) {
            db 'Database to update / create', args:1, required:true
        }
    }
}
