#!/bin/bash --login

#SBATCH --time=42:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=batch
#SBATCH -o out_ibex/%x-%j.
#SBATCH --mem=200GB

SECONDS=0


## Need to have following folders
#mkdir  10_seurat_pigeon        2_lima           5_isoseq_umicorrection  8_isoseq_collapsetranscripts  9_sqanti_filter   10_seurat_pigeon_novel  3_isoseq_bc_umi  6_isoseq_umidedup       9_pigeon_classify  1_skera                 4_isoseq_polyA   7_pbmm2                 9_sqanti                      


# ############################################ PREPROCESSING ########################################################

## 1. Read Segmentation `Skera`
skera split m84180_240702_065841/m84180_240702_065841.hifi_reads.bam \
     adapters.fasta 1_skera/m84180_240702_065841.skera.bam

#2. Remove primers 
# Remove 5' and 3' cDNA primer
module load lima

lima --isoseq --dump-clips 1_skera/m84180_240702_065841.skera.bam /ref/10x_primers.fasta 2_lima/skera_noprimer.bam

## 3. Extract UMI and cell barcodes 
module load isoseq3
isoseq3 tag 2_lima/skera_noprimer.5p--3p.bam 3_isoseq_bc_umi/skera_noprimer_taggedBCUMI.bam --design T-12U-16B

## 4. Remove poly A
module load isoseq3

isoseq3 refine --require-polya 3_isoseq_bc_umi/skera_noprimer_taggedBCUMI.bam \
    /ref/10x_primers.fasta \
    /4_isoseq_polyA/output.5p--3p.tagged.refined.bam

## 5. Barcode correction 
module load isoseq3

isoseq3 correct --barcodes /ref/barcode_3M-february-2018-REVERSE-COMPLEMENTED.txt \
   4_isoseq_polyA/output.5p--3p.tagged.refined.bam \
   5_isoseq_umicorrection/output.5p--3p.tagged.refined.corrected.bam

## 6. UMI deduplication
module load samtools
module load isoseq3

samtools sort -t CB \
   -o 5_isoseq_umicorrection/output.5p--3p.tagged.refined.corrected.sorted.bam 5_isoseq_umicorrection/output.5p--3p.tagged.refined.corrected.bam

isoseq3 groupdedup \
   5_isoseq_umicorrection/output.5p--3p.tagged.refined.corrected.sorted.bam \
   6_isoseq_umidedup/output.5p--3p.tagged.refined.corrected.sorted.dedup.bam 



7. Mapping

conda activate long_reads_prepross


pbmm2 align --preset ISOSEQ --sort 6_isoseq_umidedup/output.5p--3p.tagged.refined.corrected.sorted.dedup.bam \
   /ref/Mouse_mm39_Gencode_vM28/mouse_GRCm39.fasta \
   mapped.bam

# 8. Collapse into unique transcripts
conda activate long_reads_prepross

isoseq3 collapse ../7_pbmm2/mapped.bam collapsed.gff


# 9.1 PIGEON 
#https://isoseq.how/classification/workflow.html

#Download files we will use https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Mouse_mm39_Gencode_vM28/ tried preparong by myself but did not work :/ 

# 1 Prepare collapsed.gff
conda activate long_reads_prepross

pigeon prepare collapsed.gff

# ## 2 classify 
conda activate long_reads_prepross
pigeon classify ../8_isoseq_collapsetranscripts/collapsed.sorted.gff /ref/Mouse_mm39_Gencode_vM28/gencode.vM28.annotation.sorted.gtf \
    /ref/Mouse_mm39_Gencode_vM28/mouse_GRCm39.fasta --cage-peak /ref/Mouse_mm39_Gencode_vM28/refTSS_v3.3_mouse_coordinate.mm10.sorted.bed \
    --poly-a /ref/Mouse_mm39_Gencode_vM28/polyA.list.txt --fl ../8_isoseq_collapsetranscripts/collapsed.abundance.txt

# 3 filter artifacts default
conda activate long_reads_prepross
pigeon filter collapsed_classification.txt --isoforms ../8_isoseq_collapsetranscripts/collapsed.sorted.gff 


## Gene Saturation
conda activate long_reads_prepross

pigeon report collapsed_classification.filtered_lite_classification.txt saturation_report.txt


## Seurat Compatible File
conda activate long_reads_prepross

pigeon make-seurat --dedup /6_isoseq_umidedup/output.5p--3p.tagged.refined.corrected.sorted.dedup.fasta --group ../8_isoseq_collapsetranscripts/collapsed.group.txt \
   -d /10_seurat_pigeon \
   collapsed_classification.filtered_lite_classification.txt


##9. Annotate each transcriptt with SQANTI
conda activate SQANTI3.env

# export PYTHONPATH=$PYTHONPATH:/ibex/scratch/projects/c2169/Cecilia/isoforms/data/processed/long_reads/cDNA_Cupcake/sequence
# export PYTHONPATH=$PYTHONPATH:/ibex/scratch/projects/c2169/Cecilia/isoforms/data/processed/long_reads/cDNA_Cupcake


python SQANTI3/sqanti3_qc.py \
     ../8_isoseq_collapsetranscripts/collapsed.sorted.gff \
    /ref/Mouse_mm39_Gencode_vM28/gencode.vM28.annotation.sorted.gtf \
    /ref/Mouse_mm39_Gencode_vM28/mouse_GRCm39.fasta \
    -fl ../8_isoseq_collapsetranscripts/collapsed.abundance.txt -t 32 --report both -o annotated9 \
    --genename \
    --isoAnnotLite \



# 10. Filter artifacts
conda activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:/ibex/scratch/projects/c2169/Cecilia/isoforms/data/processed/long_reads/cDNA_Cupcake/sequence
export PYTHONPATH=$PYTHONPATH:/ibex/scratch/projects/c2169/Cecilia/isoforms/data/processed/long_reads/cDNA_Cupcake


python SQANTI3/sqanti3_filter.py rules \
     ../9_sqanti/annotated9_classification.txt \
     --isoAnnotGFF3 ../9_sqanti/annotated9.gff3 \
     --isoforms ../9_sqanti/annotated9_corrected.fasta \
     --gtf ../9_sqanti/annotated9_corrected.gtf \
     -o filtered10 \
     -d 9_sqanti_filter/ \

### Filter artifacts
grep -v 'Artifact' filtered10_RulesFilter_result_classification.txt > filtered_classification.txt




# Annotate NEED TO PREPROCESS SHORT!!!!
python cDNA_Cupcake/singlecell/link_molecule_to_celltype.py \
    /ref/bc_info.csv \
    10_seurat_pigeon/collapsed.annotated.info.csv \
    10_seurat_pigeon/dedup.annotated.celltype.csv








