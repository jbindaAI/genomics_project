#!/usr/bin/env bash
DIR=""
MSA="${DIR}/l9/msa"
OUT="${DIR}/l9/d_nieukorzenione"

N=7     # number of parallel jobs

parallel -j $N raxml-ng --threads 10
