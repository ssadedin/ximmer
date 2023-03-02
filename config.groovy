BASE="/group/bioi1/shared/tools/ximmer"
HGFA="/group/bioi1/shared/genomes/hg19/decoy/ucsc.hg19.with_decoy.fasta"

// TODO: set this to reference file containing CNVs from DGV
// You can download the required file from the UCSC DGV annotation file, "dgvMerged.txt"
// which was located at this URL on 5th feb 2015:
// http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/dgvMerged.txt.gz
// Make sure you unzip it.
DGV_CNVS="/group/bioi1/shared/tools/ximmer/eval/pipeline/tools/dgv/dgvMerged.txt.gz"

// ======================= Advanced properties =======================

TOOLS="$BASE/../../tools"

concurrency=32

// Needed by CRAM support inside anything using Picard
System.properties.reference=HGFA

GATK="$TOOLS/gatk/2.3.9"

EXCAVATOR="tools/excavator/2.0"
XHMM="$TOOLS/xhmm/trunk/build/execs"  

// These are found in the path by default, but you can 
// set specific paths here if you wish to resolve to a particular
// samtools version
SAMTOOLS="samtools"

CRAMTOOLS="$TOOLS/cramtools/2.1/cramtools.jar"

CONIFER="$TOOLS/conifer/0.2.2"

GNGS_JAR="$TOOLS/groovy-ngs-utils/current/groovy-ngs-utils.jar"

// ======================= Languages =======================

// Groovy distribution for running in the pipeline (should be the same
// as that which the bundled bpipe uses, as well as the one that
// the groovy-ngs-utils is compiled with
GROOVY="$TOOLS/groovy/2.4.6/bin/groovy"

JAVA="java"

PYTHON="/group/bioi1/shared/tools/ximmer/eval/pipeline/tools/python/miniconda/bin/python"

PERL5="perl"

