#!/bin/bash -l
#$ -l h_rt=2:00:00
#$ -N stitchpeaks
#$ -m e
#$ -o 01_stitchpeaks2.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 2

module load python2/2.7.16 
module load R
cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun/stitchpeaks/rose_strict_share/

python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/1-SS-H3K27me3/1SS27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_unique.sorted.dedup.bam

python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/2-C-H3K27me3/2C27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-C-H3K27me3/2-HSC3-C-H3K27me3_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-C-H3K27me3/2-HSC3-C-H3K27me3_unique.sorted.dedup.bam


python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/2-SS-H3K27me3/2SS27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_unique.sorted.dedup.bam

python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/CH3K27me3/C27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Cmplt-H3K27me3/Cmplt-H3K27me3_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K27me3/Cmplt-H3K27me3_unique.sorted.dedup.bam

python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/SSH3K27me3/SS27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/SS-H3K27me3/SS-H3K27me3_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K27me3/SS-H3K27me3_unique.sorted.dedup.bam


python my_rose.py -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/1-C-H3K27me3/1C27gap -w 4000 /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-C-H3K27me3.rep1/1-HSC3-C-H3K27me3.rep1_peaks.gappedPeak hg19 /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-C-H3K27me3/1-HSC3-C-H3K27me3_unique.sorted.dedup.bam
