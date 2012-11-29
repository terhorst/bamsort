#!/bin/bash -x

function bench {
  for n in 1 2 4 8 16; do
      time -o bench/$1_$n ./bamsort /data/gatk_bundle/NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam ./out.bam -n $n
  done 
}

bench python
. ~/pypy-2.2/bin/activate
bench pypy


