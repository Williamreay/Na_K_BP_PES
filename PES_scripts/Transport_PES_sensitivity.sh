#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=50:00:00

## Transport sensitivity PES

while read FILE OUTPUT;
  do
    /home/wrr239/c3202333/plink2 --pfile /home/wrr239/c3202333/Variant_QC_FINAL/Biallelic_common_merged_plink2_binary/RARE_AND_COMMON_BIALLELIC \
    --score /home/wrr239/c3202333/Na_K_UKBB_Erin_Clare/UKBB_PGS/Transport_PES_sensitivity_scores/${FILE} header-read ignore-dup-ids --out /home/wrr239/c3202333/Na_K_UKBB_Erin_Clare/UKBB_PGS/Transport_PES_sensitivity_scores/${OUTPUT};
  done < /home/wrr239/c3202333/Na_K_UKBB_Erin_Clare/UKBB_PGS/Transport_input.txt
