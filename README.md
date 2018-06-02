#!/usr/local/bin/bash

# CoNVex is an R package, developed with some unconventional methods for its use in large-scale linux clusters with large number of patient samples (genomic data)

# Execute following commands in Unix to idownload, build and install CoNVex as an R package 

# Downloadi or clone this repo:
git clone https://github.com/vijayaparthiban/convex.git

echo 'Building CoNVex package'
R CMD build CoNVex
wait
hi=`ls C*.tar.gz | sort | tail -n1`

echo Installing CoNVex package: R CMD INSTALL $hi
R CMD INSTALL $hi
wait
