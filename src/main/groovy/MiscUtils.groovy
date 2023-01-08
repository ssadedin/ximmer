import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler
import java.util.logging.Level;
import java.util.regex.Pattern
import java.nio.file.FileSystems
import java.nio.file.PathMatcher
import java.text.DateFormat;
import java.text.SimpleDateFormat
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler
import java.util.logging.Level
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import java.util.regex.Pattern
import java.util.logging.Formatter;


import groovy.util.logging.Log;;

/**
 * The default Java log former uses a format that is too verbose, so
 * replace it with something more compact.
 */
public class LogFormatter extends Formatter {
    
    private static final String lineSep = System.getProperty("line.separator");
    
    DateFormat format = new SimpleDateFormat("h:mm:ss");
    
    /**
     * A Custom format implementation that is designed for brevity.
     */
    public String format(LogRecord record) {
        
    
        String loggerName = record.getLoggerName();
        if(loggerName == null) {
            loggerName = "root";
        }
        StringBuilder output = new StringBuilder()
            .append(loggerName)
            .append("\t[")
            .append(record.threadID)
            .append("]\t")
            .append(record.getLevel()).append("\t|")
            .append(format.format(new Date(record.getMillis())))
            .append(' ')
            .append(record.getMessage()).append(' ')
            .append(lineSep);
            
        if(record.getThrown()!=null) {
            StringWriter w = new StringWriter()
            record.getThrown().printStackTrace(new PrintWriter(w))
            output.append("Exception:\n" + w.toString())
        }    
            
        return output.toString();
    }
 
}


@Log
class MiscUtils {
    
   	// See https://gist.github.com/310321
    // Licensed under the "Apache License, Version 2.0" (c) 2010
    /**
     * Returns filenames found by expanding the passed pattern which is String or
     * a List of patterns.
     * NOTE: that this pattern is not a regexp (it's closer to a shell glob).
     * NOTE: that case sensitivity depends on your system.
     *
     * <code>*</code>      Matches any file. Can be restricted by other values in
     *                     the glob pattern (same as <code>.*</code> in regexp).
     *                     <code>*</code> will match all files,
     *                     <code>c*</code> will match all files beginning with c,
     *                     <code>*c</code> will match all files ending with c.
     *
     * <code>**</code>     Matches directories recursively.
     *
     * <code>?</code>      Matches any one character. Equivalent to <code>.</code>
     *                     in a regular expression.
     *
     * <code>[set]</code>  Matches any one character in set. Behaves like character
     *                     sets in regex, including negation (<code>[^a-z]</code>).
     *
     * <code>{p,q}</code>  Matches either literal p or literal q. Matching literals
     *                     may be more than one character in length. More than two
     *                     literals may be specified. Same as alternation in regexp.
     *
     * NOTE: when matching special characters an escape is required, for example :
     * <code>"\\*"</code> or <code>"\\\\"</code>.
     *
     * NOTE: flags (e.g. case insensitive matching) are not supported.
     *
     * @see http://ruby-doc.org/core/classes/Dir.html
     * @see http://www.faqs.org/docs/abs/HTML/globbingref.html
     * @author Karol Bucek
     */
    static List<String> glob(rawPattern) {

        log.info "Globbing $rawPattern .... "

        String pattern = rawPattern
        if(System.properties['os.name'].toLowerCase().contains('windows')) {
            return windowsGlob(rawPattern)
        }
        
        if(pattern instanceof Pattern)
            return regexGlob(pattern)
        
        if ( pattern == null ) throw new IllegalArgumentException('null pattern')
        if ( pattern instanceof Collection
            || pattern instanceof Object[] ) {
            if ( pattern.size() == 0 ) return []
            return pattern.toList().sum({ glob(it) })
        }
        def base = '', path = pattern.tokenize('/')
        int i = -1, s = path.size()
        while ( ++i < s - 1 ) {
            // STOP on 'wild'-cards :
            // 1. * (equivalent to /.*/x in regexp)
            // 2. ? (equivalent to /.{1}/ in regexp)
            // 3. [set]
            // 4. {p,q}
            if ( path[i] ==~ /.*[^\\]?[\*|\?|\[|\]|\{|\}].*/ ) break
        }
        base = path[0..<i].join('/'); 
        if(pattern.startsWith("/"))
            base = "/" + base
        pattern = path[i..<s].join('/') 
        
        // a char loop over the pattern - instead of a bunch of replace() calls :
        char c; boolean curling = false; // (c) Vancouver 2010 :)
        final Closure notEscaped = { j -> // todo handling 2 escapes is enought !
            if ( j == 0 || pattern.charAt(j-1) != '\\' ) return true
            return ( j > 1 && pattern.charAt(j-2) == '\\') // [j-1] was '\\'
        }
        StringBuilder pb = new StringBuilder()
        for (i=0; i<(s = pattern.length()); i++) {
            switch (c = pattern.charAt(i)) {
                case ['.', '$'] as char[] : // escape special chars
                    pb.append('\\').append(c)
                    break
                case '?' as char : // 2. ?
                    if ( notEscaped(i) ) pb.append('.')
                    else pb.append(c)
                    break
                case '*' as char : // 1. * (or **)
                    if ( notEscaped(i) ) {
                        if ( i==s-1 || pattern.charAt(i+1) != '*' ) pb.append('.*?')
                        else (pb.append('**') && i++) // skip next *
                    }
                    else pb.append(c)
                    break
                case '{' as char : // 4. {a,bc} -> (a|bc)
                    if ( notEscaped(i) ) { pb.append('('); curling = true }
                    else pb.append(c)
                    break
                case ',' as char : // 4. {a,bc} -> (a|bc)
                    if ( notEscaped(i) && curling ) pb.append('|')
                    else pb.append(c)
                    break
                case '}' as char : // 4. {a,bc} -> (a|bc)
                    if ( notEscaped(i) && curling ) { pb.append(')'); curling = false }
                    else pb.append(c)
                    break
                default : pb.append(c)
            }
        }
        // if the last char is not a wildcard match the end :
        if ( c != '?' && c != ')' && c != ']' ) pb.append('$')
        pattern = pb.toString()
        // meh - a nice one :
        // new File('').exists() != new File(new File('').absolutePath).exists()
        final File baseFile = new File(base).getAbsoluteFile() // base might be ''
        final List fnames = [] // the result - file names
        //println "base: $base pattern: $pattern"
        if ( baseFile.exists() ) { // do not throw a FileNotFoundException
            final List matchedDirs = [ baseFile ]
            if ( (path = pattern.tokenize('/')).size() > 1 ) {
                // list and yield all dirs of the given dir :
                final Closure listDirs = { dir, yield ->
                    for ( File file : dir.listFiles() )
                        if ( file.isDirectory() ) yield.call(file, yield)
                }
                path[0..-2].each { subPattern ->
                    final boolean global = (subPattern == '**')
                    // match the dir, second param is the closure itself :
                    final Closure matchDir = { dir, self ->
                        if ( global || dir.name ==~ subPattern ) {
                            matchedDirs.add(dir)
                        }
                        if ( global ) listDirs(dir, self) // recurse
                    }
                    File[] mdirs = matchedDirs.toArray(); matchedDirs.clear()
                    for ( File mdir : mdirs ) {
                        if ( global ) matchedDirs.add(mdir)
                        listDirs( mdir, matchDir )
                    }
                }
            }
            // we used the absolute path - thus might need to remove the 'prefix' :
            s = base ? baseFile.path.lastIndexOf(base) : (baseFile.path.length() + 1)
            // add the files matching in a given directory to the result :
            final Closure addMatchingFiles = { dir, p ->
                dir.list({ pdir, name ->
                    if ( name ==~ p ) fnames << "${pdir.path}/$name".substring(s)
                    return false // we do not care about the list() return value
                } as FilenameFilter)
            }
            for (i = 0; i<matchedDirs.size(); i++) {
                // we only need the match agains the last "path"
                // aka the pattern was tokenized with '/' :
                addMatchingFiles(matchedDirs[i], path[-1])
            }
        }
        return fnames.sort()
    }

    static List<String> regexGlob(Pattern globPattern) {
        File f = new File(globPattern.toString())
        File dir = f.parentFile != null ? f.parentFile : new File(".")
        Pattern pattern = Pattern.compile(f.name)
        def result = dir.listFiles().grep { pattern.matcher(it.name).matches() }*.path        
        return result
    }

<<<<<<< HEAD
    static List<String> windowsGlob(String rawPattern) {
        File parent = new File('.')
        File patternFile = new File(rawPattern)
        File start = parent
        String pattern = rawPattern

        if(patternFile.absolute) {
            start = patternFile.absoluteFile.parentFile
            pattern = patternFile.name
            return new FileNameFinder().getFileNames(start.absolutePath, pattern)
        }

        PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:" + pattern);
        List results = []
        start.eachFileMatch(matcher) { 
            results.add(it.path)
        }
        return results 
    }

    /**
     * Set up simple logging to work in a sane way
     *
     * @param path
     * @return
     */
    public static configureSimpleLogging(String path) {
        Logger parentLog = log.parent
        parentLog.getHandlers().each { parentLog.removeHandler(it) }
        FileHandler fh = new FileHandler(path)
        fh.setFormatter(new LogFormatter())
        parentLog.addHandler(fh)
    }

    public static void configureVerboseLogging() {
        ConsoleHandler console = new ConsoleHandler()
        console.setFormatter(new LogFormatter())
        console.setLevel(Level.FINE)
        log.getParent().addHandler(console)
    }

    final static int TEMP_DIR_ATTEMPTS = 10

    /**
     * Create a temporary directory 
     * (Based on Guava library)
     */
    public static File createTempDir() {
        File baseDir = new File(System.getProperty("java.io.tmpdir"));
        String baseName = Long.toString(System.nanoTime()) + "-";

        for (int counter = 0; counter < TEMP_DIR_ATTEMPTS; counter++) {
            File tempDir = new File(baseDir, baseName + counter);
            if (tempDir.mkdir()) {
                return tempDir;
            }
        }
        throw new IllegalStateException("Failed to create directory within $TEMP_DIR_ATTEMPTS")
    }
}
