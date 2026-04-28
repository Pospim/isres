#!/bin/bash

findMotifs.pl raw/hs_ifn1_isg_gal_orth_2kb_up.fa fasta homer_hs_ifn1_ort_2kb_cpg \
  -fastaBg raw/bcg/chicken_bcg_promotors_20260320.fa \
  -len 8,10,12,14 \
  -p 12 \
  -S 16 \
  -cpg
