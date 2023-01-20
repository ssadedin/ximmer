BASE="."
TOOLS="$BASE/tools"
PADDING=25

width = 100


hd = { msg ->
	println("\n"+(" " + msg + " ").center(width,"=") + "\n")
}

file = { new File(it) }

stdin = System.in.newReader()

def ask(String msg, Closure c) {
    print msg + "? (y/n) "
    answer = stdin.readLine()
    if(answer == "y") {
        c()
    }
}

def sh(String script, Map env = [:]) {
    
    String dir = "."
    
    String rTempDir = MiscUtils.createTempDir().absolutePath
    File tempScript = new File(rTempDir, "XimmerInstallScript.sh")
    tempScript.setExecutable(true)
    tempScript.text = script 
    ProcessBuilder pb = new ProcessBuilder([
       "bash", tempScript.absolutePath
    ] as String[]).directory(new File(dir))
        
    Map pbEnv = pb.environment()
    env.each { k,v ->
        pbEnv[k] = String.valueOf(v)
    }
                      
    Process process = pb.start()
    process.consumeProcessOutput(System.out, System.err)
            
    return process.waitFor()
}

ant = new AntBuilder()

class ClosureScript extends Script {
    Closure closure
    def run() {
        closure.resolveStrategy = Closure.DELEGATE_FIRST
        closure.delegate = this
        closure.call()
    }
}

tools = {
    xhmm {
        
        dir = "$BASE/pipeline/tools/xhmm"
        
        exists = {
            file("$BASE/pipeline/tools/xhmm/trunk/xhmm").exists()
        }
        
        doinstall = """
            mkdir -p $BASE/pipeline/tools/xhmm
            cd $BASE/pipeline/tools/xhmm
            curl -o master.zip https://bitbucket.org/statgen/xhmm/get/master.zip
            unzip master.zip
            mv statgen* trunk
            cd trunk
            echo
            echo "Patching XHMM to fix compile errors on newer C++ ..."
            echo
            sed -i.bak  's/noexcept.*{/{/g'   ./sources/hmm++/sources/include/LAPACKvector.hpp  ./sources/hmm++/sources/include/VectorOnDisk.hpp 
            make
        """
        
        check { 
            is_xhmm_executable = "$BASE/pipeline/tools/xhmm/trunk/xhmm --version | grep -q 'xhmm [0-9]'"
        }
        
        license = """
            XHMM licensing is unclear; you should review the statements on the 
            XHMM web site [https://atgu.mgh.harvard.edu/xhmm/] and on the user group
            forum https://groups.google.com/a/broadinstitute.org/forum/#!topic/xhmm-users/zNqbOKFwtUs
            to decide if you are eligible to use XHMM.
        """
    }
}

ConfigSlurper slurper = new ConfigSlurper()
slurper.setBinding(this.binding.variables)
ConfigObject config = slurper.parse(new ClosureScript(closure:tools))


println "*" * width
println "Ximmer Install Script"
println "*" * width
println ""
println """
This script helps to install Ximmer and check that all the necessary tools are 
available and configured correctly.
"""

errors = []

def install_tool(String toolName, ConfigObject tool) {
    
    hd toolName
    if(!("check" in tool))  {
        println "Nothing to do for $toolName"
    }
    
    if(!tool.exists()) {
        
        println "$toolName is not found / compiled / installed."
        
        if("doinstall" in tool) {
            
            ask("Do you want to install $toolName") {
            
                if("license" in tool) {
                    println tool.license
                    println ""
                    ask("Do you wish to continue downloading / installing $toolName") {
                        sh tool.doinstall
                    }
                }
                else {
                    sh tool.doinstall
                }
            }
        }
        else {
            error "Tool $toolName is not installed, and Ximmer does not know how to install it. Please install $toolName manually."
            errors << error
            return
        }
    }
    
    for(Map.Entry check in tool.check) {
        
        boolean result = false
        if(check.value instanceof String || check.value instanceof GString) {
            try {
                result = (sh(check.value) == 0)
            }
            catch(Exception e) {
                // failed to execute command
                e.printStackTrace()
            }
        }
        else
        if(tool.check.value instanceof Closure){
            result = tool.check.value()
        }
        else {
            System.err.println "ERROR: install script error - was expecting $check.key to refer to either a String or a Closure but it was a ${check.value.class.name}"
            System.exit(1)
        }
            
        println check.key.padRight(PADDING) + ": " + (result?"yes":"no")
        
        if(!result) {
            def error = "ERROR: check $check.key for tool $toolName failed"
            errors << error
        }
    }
}

for(tool in config) {
    install_tool(tool.key, tool.value)    
}

hd "Finished"

if(!errors) {
    println "Success - everything checks out!"
}
else {
    println "Some problems occurred: \n"
    
    errors.each {
        println "* " + it
    }
    
    println ""
    println "Please resolve these before trying to run Ximmer!"
}


