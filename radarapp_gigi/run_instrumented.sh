#!/bin/bash
source /apps/riscv/extrae/3.8.3/etc/extrae.sh
export LD_LIBRARY_PATH=/apps/riscv/papi-like/development/lib:$LD_LIBRARY_PATH
export EXTRAE_CONFIG_FILE=extrae.xml
LD_PRELOAD=/apps/riscv/extrae/3.8.3/lib/libseqtrace-3.8.3.so "$@"

