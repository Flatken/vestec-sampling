#!/bin/bash

module() { eval `/usr/bin/modulecmd bash $*`; }

export MODULEPATH=/tools/modulesystem/modulefiles
export LC_ALL=C

module purge
module load clang/clang-4.0.0

SRC_DIR="$( cd "$( dirname "$0" )" && pwd )"

FILES=`find $SRC_DIR | grep ".cpp\|.hpp\|.inl"`

for FILE in $FILES; do
    clang-format -i $FILE
done
