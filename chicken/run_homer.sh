#!/bin/bash

findMotifs.pl raw/hs_orth_isgs_ensmbl_list.txt chicken homer_hs_orth_builtinbg \
  -start -2000 -end 0 \
  -len 8,10,12 \
  -p 12 \
