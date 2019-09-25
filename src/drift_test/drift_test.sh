#!/bin/sh

HTSDIR=../htslib
BCFDIR=../bcftools
export LD_LIBRARY_PATH=$HTSDIR:$LD_LIBRARY_PATH
export BCFTOOLS_PLUGINS=./
#valgrind --tool=exp-sgcheck
#valgrind --leak-check=full --show-leak-kinds=all \
$BCFDIR/bcftools +af_drift_test $@
