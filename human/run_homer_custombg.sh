#!/bin/bash

findMotifs.pl raw/ifn1_2kb_up.fa fasta homer_2kb_custombg \
  -fastaBg raw/bcg/hs_2kb_bcg_20260328.fa \
  -len 8,10,12,14 \
  -p 12 \
  -S 16 \
