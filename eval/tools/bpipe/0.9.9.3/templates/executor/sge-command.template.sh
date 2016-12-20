#!/bin/sh
#\$ -wd ${bpipe.Runner.canonicalRunDirectory}
#\$ -N "$name"
#\$ -terse
#\$ -o $jobDir/$CMD_OUT_FILENAME
#\$ -e $jobDir/$CMD_ERR_FILENAME
#\$ -notify
<%if(config?.walltime) {%>#\$ -l h_rt=${config.walltime} <%}
procs=null
if(config.containsKey('procs') && config.procs.toString().isInteger()) {
    procs = config.procs.toString().toInteger()
}
if(procs>1) {
    if(config.containsKey('sge_pe')) {%>
#\$ -pe ${config.sge_pe} ${config.procs} <%
    }
    else {%>
#\$ -l slots=${config.procs} <%
    }
}
if(config?.memory) {
    mem_param = 'virtual_free'
    if(config.containsKey('mem_param')) {
        mem_param = config.mem_param 
    } 
    if(procs>1 && config.containsKey('smp_mem_param' )) {
        mem_param = config.smp_mem_param 
    }
}%>
#\$ -l ${mem_param}=${config.memory}

# Executed usinq queue ${config.queue}

function signal_handler() {
    echo 255 > $jobDir/$CMD_EXIT_FILENAME;
}

export signal_handler

trap signal_handler HUP INT TERM QUIT KILL STOP USR1 USR2

<%if(config.containsKey('use')) {%>
<%if(config.containsKey('use_prefix')) {%> ${config.use_prefix} <%}%>
use ${config.use}
<%}%>

(
bash -e ${CMD_FILENAME.absolutePath}
) &

wait \$!

result=\$?
echo -n \$result > $jobDir/$CMD_EXIT_FILENAME
exit \$result
