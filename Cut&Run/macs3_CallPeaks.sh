#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N callpeaksbroad
#$ -m e
#$ -o 09_callpeaksbroad.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 16

cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun
module load python3/3.8.6
module load macs/3.0.0a6
#IgGs (for use as blacklist in heatmaps only)
macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Ctrl_IgG/ -n Ctrl_IgG

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-IgG/ -n 1-IgG

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-IgG/ -n 2-IgG

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K27me3/Cmplt-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Cmplt-H3K27me3/ -n Cmplt-H3K27me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Cmplt-H3K4me3/Cmplt-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Cmplt-H3K4me3/ -n Cmplt-H3K4me3   

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K27me3/SS-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/SS-H3K27me3/ -n SS-H3K27me3   

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/SS-H3K4me3/SS-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/20240215_SJ_SerineStarvationHSC3/20240723_normalized/align/Ctrl-IgG/Ctrl-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/SS-H3K4me3/ -n SS-H3K4me3  


macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-C-H3K27me3/1-HSC3-C-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-C-H3K27me3.rep1/ -n 1-HSC3-C-H3K27me3.rep1

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-SS-H3K4me3/1-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam  --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K4me3/ -n 1-HSC3-SS-H3K4me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-C-H3K4me3/2-HSC3-C-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-C-H3K4me3/ -n 2-HSC3-C-H3K4me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K27me3/ -n 2-HSC3-SS-H3K27me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-C-H3K4me3/1-HSC3-C-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-C-H3K4me3/ -n 1-HSC3-C-H3K4me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/1-HSC3-IgG/1-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K27me3/ -n 1-HSC3-SS-H3K27me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-C-H3K27me3/2-HSC3-C-H3K27me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-C-H3K27me3/ -n 2-HSC3-C-H3K27me3

macs3 callpeak -B -q 0.05 --broad --keep-dup 1 -g hs -f BAMPE --cutoff-analysis --extsize 146 --nomodel -t /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-SS-H3K4me3/2-HSC3-SS-H3K4me3_unique.sorted.dedup.bam -c /restricted/projectnb/montilab-p/CBMrepositoryData/otherStudies/CutNRun/240801_SJ12268_chips/align/2-HSC3-IgG/2-HSC3-IgG_unique.sorted.dedup.bam --outdir /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K4me3/ -n 2-HSC3-SS-H3K4me3
