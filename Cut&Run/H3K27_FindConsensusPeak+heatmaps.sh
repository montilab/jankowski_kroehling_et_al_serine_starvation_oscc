#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N h3k27hm
#$ -m e
#$ -o 25_h3k27heatmaps.log
#$ -j y
#$ -P gsc-p
#$ -pe omp 4

cd /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/scripts/cutnrun
module load python3/3.8.6
module load deeptools/3.5.1
module load bedtools

##first for ss
##first get merged peaks
cat /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_peaks.broadPeak | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 -c 4,5 -o distinct,max > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27_merged_peaks.bed 

#get peaks found in multiple files
bedtools multiinter -i /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-HSC3-SS-H3K27me3/1-HSC3-SS-H3K27me3_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-HSC3-SS-H3K27me3/2-HSC3-SS-H3K27me3_peaks.broadPeak | sort -k1,1 -k2,2n > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27peak_counts.bed

#get peaks in at least 2 files
awk '($4>=2) {print $0}' /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27peak_counts.bed > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27peak_counts_in2.bed

#pull peaks in multiple files from merged peak file
bedtools intersect -a /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27_merged_peaks.bed -b /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27peak_counts_in2.bed -wa -u -sorted > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27_consensusPeaks.bed


###second for control
########control H3K27
cat /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/CH3K27me3/C27gap_Gateway_Enhancers.bed /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/1-C-H3K27me3/1C27gap_Gateway_Enhancers.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 -c 4,5 -o distinct,max > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_merged_peaks.bed

#get peaks found in multiple files
bedtools multiinter -i /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/CH3K27me3/C27gap_Gateway_Enhancers.bed /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/stitchpeaks/1-C-H3K27me3/1C27gap_Gateway_Enhancers.bed | sort -k1,1 -k2,2n > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitchpeak_counts.bed

#get peaks in at least 2 files
awk '($4>=2) {print $0}' /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitchpeak_counts.bed > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitchpeak_counts.bed_in2.bed

#pull peaks in multiple files from merged peak file
bedtools intersect -a /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_merged_peaks.bed -b /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitchpeak_counts.bed_in2.bed -wa -u -sorted > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_consensusPeaks.bed

###third for IgGs
cat /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Ctrl_IgG/Ctrl_IgG_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-IgG/1-IgG_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-IgG/2-IgG_peaks.broadPeak | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 -c 4,5 -o distinct,max > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgG_merged_peaks.bed

#get peaks found in multiple files
bedtools multiinter -i /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/Ctrl_IgG/Ctrl_IgG_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/1-IgG/1-IgG_peaks.broadPeak /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/2-IgG/2-IgG_peaks.broadPeak | sort -k1,1 -k2,2n > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgG_peak_counts.bed

#get peaks in at least 2 files
awk '($4>=2) {print $0}' /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgG_peak_counts.bed > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgGpeak_counts.bed_in2.bed

#pull peaks in multiple files from merged peak file
bedtools intersect -a /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgG_merged_peaks.bed -b /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgGpeak_counts.bed_in2.bed -wa -u -sorted > /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/IgG_consensusPeaks.bed

######heatmaps
##Complete 
###H3K27 Complete
# Compute matrix
computeMatrix scale-regions -S /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Ctrl-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-C-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K27me3norm.bw -R /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_consensusPeaks.bed -b 5000 -a 5000 -m 5000 --maxThreshold 63000 --binSize 50 -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/1CmpltH3K27stitch_CP_matrixMT.gz

# Plot heatmap
plotHeatmap -m /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/1CmpltH3K27stitch_CP_matrixMT.gz -out /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/20241122_Cmplt_H3K27_MT.png --colorMap RdYlBu --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --zMin -3 --zMax 3 --refPointLabel "Peak Center" --regionsLabel "Cmplt H3K27 Stitched Consensus Peaks" --plotTitle "CUT&RUN Signal at Cmplt H3K27 Stitched Consensus Peaks"


###H3K27 Complete reference point
# Compute matrix #######NOT THIS ONE
#computeMatrix reference-point -S /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-C-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K27me3norm.bw -R /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/Cmplt_H3K27stitch_consensusPeaks.bed --referencePoint center -b 50000 -a 50000 --maxThreshold 63000 --binSize 50 -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/1CmpltH3K27stitch_CP_matrixREF.gz

# Plot heatmap
#plotHeatmap -m /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/1CmpltH3K27stitch_CP_matrixREF.gz -out /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/1CmpltH3K27stitch_CPheatmapREFBL.png --colorMap RdYlBu --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --zMin -3 --zMax 3 --refPointLabel "Peak Center" --regionsLabel "Cmplt H3K27 Stitched Consensus Peaks" --plotTitle "CUT&RUN Signal at Cmplt H3K27 Stitched Consensus Peaks REFPOINT"



###H3K27 SS
# Compute matrix
computeMatrix reference-point -S /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-C-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-IgGnorm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/Cmplt-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/1-HSC3-SS-H3K27me3norm.bw /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/bamCoveragenorm/2-HSC3-SS-H3K27me3norm.bw -R /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/consensuspeaks/SS_H3K27_consensusPeaks.bed --referencePoint center -b 50000 -a 50000 --maxThreshold 6500 --binSize 50 -o /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/SSH3K27broad_CP_matrixMT.gz

# Plot heatmap
plotHeatmap -m /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/SSH3K27broad_CP_matrixMT.gz -out /restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/cutnrun/peaksnoSMPR_broad/visualizations/20241122_SSH3K27broad_CPheatmap50kMT.png --colorMap RdYlBu --whatToShow 'plot, heatmap and colorbar' --heatmapHeight 15 --heatmapWidth 4 --zMin -3 --zMax 3 --refPointLabel "Peak Center" --regionsLabel "SS H3K27 Broad Consensus Peaks" --plotTitle "CUT&RUN Signal at SS H3K27 Broad Consensus Peaks"
