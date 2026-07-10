#!/bin/bash -l
#$ -l h_rt=20:00:00
#$ -N bamcov
#$ -m e
#$ -o 03_bamCoverageNormalization.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 8

#use scale factor based on westerns normalized to H3 instead of RPGC normalization, as recommended by a reviewer.

cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun/2026_01_revisions/20260409_bw_westernscale
module load python3/3.8.6
module load deeptools/3.5.1

#https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K27me3/Cmplt-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/Cmplt-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.848462547

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K4me3/Cmplt-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/Cmplt-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K27me3/SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/SS-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.380283105

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K4me3/SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/SS-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.525847916

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/1-HSC3-C-H3K27me3/1-HSC3-C-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/1-HSC3-C-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.848462547

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/1-HSC3-SS-H3K4me3/1-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/1-HSC3-SS-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.525847916

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/2-HSC3-C-H3K4me3/2-HSC3-C-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/2-HSC3-C-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/2-HSC3-SS-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.380283105

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/1-HSC3-C-H3K4me3/1-HSC3-C-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/1-HSC3-C-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/1-HSC3-SS-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.380283105

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/2-HSC3-C-H3K27me3/2-HSC3-C-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/2-HSC3-C-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.848462547

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240801_SJ12268_chips/align/2-HSC3-SS-H3K4me3/2-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/2-HSC3-SS-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.525847916

#new reads
bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/A-Ctrl-H3K4me3/A-Ctrl-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/A-Ctrl-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/A-S_G-H3K4me3/A-S_G-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/A-S_G-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.525847916

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/B-Ctrl-H3K4me3/B-Ctrl-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/B-Ctrl-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/B-S_G-H3K4me3/B-S_G-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/B-S_G-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.525847916

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/B-Ctrl-H3K27me3/B-Ctrl-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/B-Ctrl-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.848462547

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/B-S_G-H3K27me3/B-S_G-H3K27me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/B-S_G-H3K27me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.380283105

bamCoverage --bam /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20250703_SJ_CutnRun_2ndRound/processed-data/250613_processed_fromCore/PreprocessedBams/C-Ctrl-H3K4me3/C-Ctrl-H3K4me3_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/20260409_bamCovnorm_Scale/C-Ctrl-H3K4me3norm.bw --binSize 10 --normalizeUsing None --scaleFactor 0.761396212
