#!/bin/bash
# Written By Schaudge King
# 2022-07-26
# build the bcl2fastq from cmake environment (e.g. CLion IDE)

source_dir=$(dirname $0)
mkdir -p ${source_dir}/build && cd ${source_dir}/build

cmake -G "CodeBlocks - Unix Makefiles" -DBCL2FASTQ_NAME_SHORT:STRING="bcl2fastq" -DBCL2FASTQ_NAME_LONG:STRING="BCL to FASTQ file converter" -DBCL2FASTQ_COPYRIGHT:STRING="Copyright (c) 2007-2017 Illumina, Inc." \
-DBCL2FASTQ_VERSION_MAJOR:STRING="2" -DBCL2FASTQ_VERSION_MINOR:STRING="21" -DBCL2FASTQ_VERSION_PATCH:STRING="0" -DBCL2FASTQ_VERSION_BUILD:STRING="422" -DBCL2FASTQ_VERSION:STRING="2.21.0.422" \
-DBCL2FASTQ_PREFIX:PATH="/usr/local" -DBCL2FASTQ_EXEC_PREFIX:PATH="" -DBCL2FASTQ_BINDIR:PATH="" -DBCL2FASTQ_LIBDIR:PATH="" -DBCL2FASTQ_LIBEXECDIR:PATH="" \
-DBCL2FASTQ_INCLUDEDIR:PATH="" -DBCL2FASTQ_DATADIR:PATH="" -DBCL2FASTQ_DOCDIR:PATH="" -DBCL2FASTQ_MANDIR:PATH="" -DCMAKE_INSTALL_PREFIX:PATH="/usr/local" \
-DCMAKE_C_FLAGS_DEBUG="-g -O0" -DCMAKE_CXX_FLAGS_DEBUG="-g -O0" -DCMAKE_PARALLEL:STRING="1" ..

make

### copy the xml template for bcl2fastq statistics
#mkdir -p cxx/share
#cp -r ../xsl ../css cxx/share/
