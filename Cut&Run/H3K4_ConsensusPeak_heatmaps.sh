#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N h3k4hm
#$ -m e
#$ -o 26_h3k4heatmaps.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 4

cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun
module load python3/3.8.6
module load deeptools/3.5.1
module load bedtools

#########SS H3K4
cat /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K4me3/1-HSC3-SS-H3K4me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K4me3/2-HSC3-SS-H3K4me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/SS-H3K4me3/SS-H3K4me3_peaks.broadPeak | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 -c 4,5 -o distinct,max > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4_merged_peaks.bed 

#get peaks found in multiple files
bedtools multiinter -i /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K4me3/1-HSC3-SS-H3K4me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K4me3/2-HSC3-SS-H3K4me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/SS-H3K4me3/SS-H3K4me3_peaks.broadPeak > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4peak_counts.bed


#get peaks in at least 2 files
awk '($4>=2) {print $0}' /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4peak_counts.bed > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4peak_counts_in2.bed

#pull peaks in multiple files from merged peak file
bedtools intersect -a /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4_merged_peaks.bed -b /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4peak_counts_in2.bed -wa -u -sorted > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4_consensusPeaks.bed



#########heatmaps
#SS consensus H3K4
# Compute matrix
computeMatrix reference-point -S /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/SS-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K4me3norm.bw -R /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K4_consensusPeaks.bed --referencePoint center -b 3000 -a 3000 --binSize 50 -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/SS_Consensus_H3K4_CP_allsamples_matrix.gz


# Plot heatmap
plotHeatmap -m /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/SS_Consensus_H3K4_CP_allsamples_matrix.gz -out /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/20241122_SSconH3K4_CP_heatmap.png --colorMap RdYlBu --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --zMin -3 --zMax 3 --refPointLabel "Peak Center" --regionsLabel "SS consensus H3K4 Consensus Peaks (>=2reps)" --plotTitle "CUT&RUN Signal at SS H3K4 Consensus Peaks"


#Complete broad peaks H3K4
# Compute matrix
computeMatrix reference-point -S /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/SS-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K4me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K4me3norm.bw -R /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Cmplt-H3K4me3/Cmplt-H3K4me3_peaks.broadPeak --referencePoint center -b 3000 -a 3000 --binSize 50 -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/Complete_H3K4_CP_allsamples_matrix.gz

# Plot heatmap
plotHeatmap -m /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/Complete_H3K4_CP_allsamples_matrix.gz -out /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/CompleteH3K4_CP_heatmap.png --colorMap RdYlBu --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --zMin -3 --zMax 3 --refPointLabel "Peak Center" --regionsLabel "C consensus H3K4 Peaks (from 1 sample)" --plotTitle "CUT&RUN Signal at Complete H3K4 Peaks"

