# Ximmer Configuration

## Introduction

Ximmer is designed to run "out of the box" with some reasonably sensible defaults.
Still, there is some unavoidable configuration if you want it to work well for you.

## Set paths for Python

Ximmer requires Python 2.7 and depends on a specific set of libraries.
If these libraries are already available, Ximmer will find them and use them
automatically. If not, it will try to install them using `pip` or `conda`. If
Python 2.7 is not your default Python version, edit the following two files
to point Ximmer to where it can find a valid Python 2.7 installation. This can
be a VirtualEnv environment, which can help to allow Ximmer to install the 
libraries it needs without interfering with other applications.

To set a custom Python executable, set the path in these two files:

 * &lt;ximmer install&gt;/eval/pipeline/config.groovy
 * &lt;ximmer install&gt;/eval/pipeline/bpipe.config

## Set up Computational Environment

By default, Ximmer runs all its operations on the same computer that you run Ximmer 
on. Under the hood, however, Ximmer uses [Bpipe](http://bpipe.org) which can interface
with computational clusters to run jobs. If you would like to run Ximmer in such a 
way, edit <ximmer install>/eval/pipeline/bpipe.config and adjust the settings in there
accordingly. See the [Bpipe Documentation](http://docs.bpipe.org/Guides/ResourceManagers/) 
for information about the right format and settings to use.

If your environment uses modules or dotkit to control versions of tools, you can also
set those in bpipe.config, or alternatively, you can change what Ximmer will use
directly by editing eval/pipeline/config.groovy.

## Set Concurrency to Use

By default Ximmer will try to use only 2 cores on the computer it runs on. If you have 
more cores available (particularly, if you are using a computational cluster), add 
a line such as the example below to your simulation/analysis config.groovy:

```
concurrency=32
```

The above will allow Ximmer to use a total of 32 cores altogether. Note these may be 
distributed across multiple commands. If using a computational cluster, set it to the
number of cores available in the whole cluster, or about 120, whichever is smaller.

## Configuring an Analysis / Simulation

See the [Analyses](analyses.md) section for how to configure an individual simulation 
or anlaysis run.
