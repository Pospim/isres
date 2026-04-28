#!/bin/bash

findMotifs.pl raw/ifn1.ensg.list human homer_ifn1_2k_builtinbg_cpg \
  -start -2000 -end 0 \
  -len 8,10,12 \
  -p 12 \
  -cpg
