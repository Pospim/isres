#!/bin/bash

findMotifs.pl chicken_2000_upstream.fa fasta motif_hits \
  -find homerMotifs.all.motifs \
  > motif_hits/homerMotifs.all_hits.txt
