#!/bin/bash

findMotifs.pl raw/hs_ifn1_upreg_duck_orth_2kb_up.fa fasta homer_hs_ifn1_orth_2k_custombg \
  -fastaBg raw/bcg/raw_2kb_bcg_20260328.fa \
  -len 10,12,14 \
  -p 12 \
  -S 16
