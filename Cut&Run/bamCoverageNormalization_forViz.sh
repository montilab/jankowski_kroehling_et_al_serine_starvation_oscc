#!/bin/bash -l
#$ -l h_rt=20:00:00
#$ -N bamcov
#$ -m e
#$ -o 03_bamCoverageNormalization.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 8

#create normalized bw files from bams using bamCoverage, this is an alternative to the SMPR normalization that can be used during macs3 callpeaks. I want to compare the two normalized files.

cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun
module load python3/3.8.6
module load deeptools/3.5.1

#https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K27me3/Cmplt-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K4me3/Cmplt-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K27me3/SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/SS-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K4me3/SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/SS-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-C-H3K27me3/1-HSC3-C-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-C-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-SS-H3K4me3/1-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-C-H3K4me3/2-HSC3-C-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-C-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-C-H3K4me3/1-HSC3-C-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-C-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-C-H3K27me3/2-HSC3-C-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-C-H3K27me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Ctrl-IgGnorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-IgGnorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-SS-H3K4me3/2-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K4me3norm.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX --extendReads
