#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ];
then
    echo
    echo "Need three arguements <assay ID> <.bed file assay full path> <CXP URL>"
    echo
    exit 1
fi

set -euo pipefail
set -x

ID="$1"
FULLPATH="$2"
XIMMER_HOME="$(realpath $(dirname $0 | xargs dirname))"

echo "$XIMMER_HOME"
CP="$XIMMER_HOME/tools/groovy-ngs-utils/1.0.9/groovy-ngs-utils.jar:$XIMMER_HOME/src/main/groovy:$XIMMER_HOME/src/main/resources:$XIMMER_HOME/src/main/js"

echo '{
    "id": null,
    "metadata": {},
    "identifier": "'$ID'",
    "alt_id": "'$ID'",
    "fullpath": "'$FULLPATH'",
    "version": 1,
    "default_reference_panel": null
}' | "$XIMMER_HOME/tools/groovy/2.5.13/bin/groovy" -cp "$CP" "$XIMMER_HOME/src/main/groovy/RawPostToCXP.groovy" \
    -url "$3/assay/"
