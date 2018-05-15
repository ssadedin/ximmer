# Installation

This section gives detailed guidance about how to get Ximmer set up and working.

There are two ways to set up Ximmer:

 * Natively - using libraries, tools and packages directly on your computer
 * Using Docker

Each approach is described below in turn.

## Installation and Requirements (Native)

To make Ximmer easier to use we have included support to automatically download
and build the tools and other dependencies it needs. You should make sure
before starting that you have at minimum the following requirements:

 * Java 1.7+ 
 * Python 2.7, preferably the Anaconda installation (you will need
   support for pandas, numpy and other computational libraries that
   are not always easy to compile in a vanilla installation).
 * A recent version of R (3.4.x as of writing)
 * 24GB of RAM
   
Ideally, these should all be directly accessible from your environment. 
If necessary, you can specify custom locations for them in the configuration 
file.

You should also make sure you have internet access while doing the installation
because Ximmer will try to download some components. It may be necessary to set 
the "http_proxy" environment variable if your network uses a proxy.


### Linux Packages

When installing on linux, you may find that the Python and R libraries that
Ximmer needs depend on the following Linux packages being installed:

 * Debian / Ubuntu : libxml2-dev libcurl4-openssl-dev libblas-dev liblapack-dev libhdf5-dev libssl-dev libmariadb-client-lgpl-dev
 * CentOS : libxml2-devel blas-dev lapack-dev libcurl-devel hdf5-devel

### Consider: Set up Python VirtualEnv Environment

The installer will attempt to install various python packages at specific
versions.  If you don't have administrative permissions to install python
packages, or if you already have conflicting versions installed, this could
create a problem.  If you have Python virtual-env installed, you can avoid
these problems by setting up a virtual-env environment to use with Ximmer. For
example:

```
virtualenv ximmer-python
source ximmer-python/bin/activate
```

### Run Installer

Ximmer includes an installer script to help set up and configure
it for basic operation. To get started:

```
git clone https://github.com/ssadedin/ximmer.git
cd ximmer
./bin/install
```


## Installation and Requirements (Docker)

If you use Docker, you can install and run Ximmer simply by building it
from the Docker file. In this case the primary requirement is that you do 
this on a computer with enough RAM, we recommend at least 24GB of RAM for 
the analysis pipeline to run successfully.

To build the Docker image use:

```bash
git clone https://github.com/ssadedin/ximmer.git
cd ximmer/docker
docker build -t ximmer . 
```

Note: if your machine is behind a proxy, you can provide it ike this:

```bash
git clone https://github.com/ssadedin/ximmer.git
cd ximmer/docker
docker build --build-arg http_proxy='http://proxy.host.com:proxy-port' -t ximmer . 
```

Although running inside docker is essentially the same as running outside,
 there are a few special considerations to take care of.  See 
 [Running inside Docker](docker.md) for tips on how to run inside Docker.

## Running Ximmer

See information about how to run an analysis using Ximmer in the [Running](running.md) 
section, and how to configure CNV simulation in [Simulations](simulations.md).
