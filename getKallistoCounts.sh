#!/bin/bash

#this script thanks to a design file, will give for each sample of each condition, the counts with kallisto.
#the results will be in the output directory that you have set
#each sample has its own sub-directory, and in each, the counts are in *_abundance.tsv


############ inputs #############
#################################

#the design should be something like this (tab-separated) :
# sample_name		condition
# A	cond1
# B	cond1
# O	cond2
# P	cond2
# ...
#design="/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/MAT1_overexp_EPIMES/test1.tsv"
#design="/home/marcgabriel/Desktop/test1.tsv"
#design="/media/marcgabriel/homeborn/Dominika_miRNA_seq/kallisto_design.tsv"
#design="/media/marcgabriel/homeborn/Dominika_snoRNA_seq/kallisto_design.tsv"
#design="/media/marcgabriel/saylar4/Nouritza_TNBClncRNA/kept_design.tsv"
#design="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/kept_design.tsv"
#design="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/kallisto_design2.tsv"
#design="/media/marcgabriel/saylar2/Nouritza_human_cell_lines/kallisto_design.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode_and_repeats/design_tumor.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode_and_repeats/design_FFPE.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/kallisto_design_transcripts.tsv"

#design="/media/marcgabriel/saylar5/prostate_cell_lines/design.txt"

#design="/media/marcgabriel/saylar5/tissues_to_use_as_control_for_urines/tissues_design.tsv"

#design="/media/marcgabriel/saylar5/nuclear_cyto_prostate_samples_SRP199212/design.tsv"
#design="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/kallisto_design_batch1And2_22_01_2021.tsv"
#design="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/kallisto_transcripts_counts_sample.tsv"
#design="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/Ciriquant_design.tsv"
#design="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/kallisto_counts_gencode32_repeats_PNT2_PC3/kallisto_design_PNT2_PC3.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/all_ORFS/TransDecoder_output/uEVs_FFPE_celllines_and_associated_EVs_design.tsv"
#design="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/all_ORFS/TransDecoder_output/kallisto_counts_on_genes/riboseq_design.tsv"
#design="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/urinary_EVs_clean_files_11_06_2021/kallisto_counts_6normaluEVs_6tumoruEVssamples/design_6NormaluEVs_6tumoruEVs.tsv"

#design="/media/marcgabriel/saylar8/plasticity_analysis/fastq_files/39_samples_analysis_kallisto_design.tsv"
#design="/media/marcgabriel/saylar8/plasticity_analysis/fastq_files/samples_analysis_kallisto_design_early_rec_vs_prolonged_surv.tsv"
#design="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/kallisto_counts_gencode32_repeats_PNT2_PC3/kallisto_updated_design.tsv"

#design="/media/marcgabriel/saylar9/circRNAs_20_uEV_samples_JEV_revision/circRNA_design_20_uEV_samples.tsv"

#design="/media/marcgabriel/saylar8/Nouritza_EV_analysis/fastq_files/kallisto_desig.tsv"


#SRR10261665_GSM4117331_RNAAtlas045_total_RNA & SRR10266661_GSM4117305_RNAAtlas036_total_RNA have an issue : all counts are 0 -> presence of adapters ? wrong lib type ?
#design="/media/marcgabriel/saylar9/Nouritza_breast_cancer_cell_lines_P_mestdagh_total_and_polyA/kallisto_design.tsv"

#design="/media/marcgabriel/saylar9/Nouritza_human_TNBClncRNA_validation_batchPart2/kallisto_design_batch_1_2_3.tsv"

#design="/media/marcgabriel/saylar11/Nouritza_TNBCdocPDX/fastq_files/kallisto_design.tsv"

#design="/media/marcgabriel/saylar10/Nouritza_TNBC_cell_lines/featurecounts_design.tsv"

#design="/media/marcgabriel/saylar14/GTex_missing_DG_samples/449_missing_DG_samples/611_and_missing_DG_and_upto50_desc_unique_no4samplesWithBadRin.tsv"

#design="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/kallisto_transcripts_counts_sample.tsv"

#design="/media/marcgabriel/saylar11/Rocco_RIPseq_analysis/fastq_files/kallisto_design.tsv"


#design="/media/marcgabriel/saylar9/urines_files_new_batches_04_04_2023/kallisto_design.tsv"

design="/media/marcgabriel/saylar5/LNCaP_data/Rocco_LNCaP_GTex_kallisto_design.tsv"

#reads path
#reads_path="/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/MAT1_overexp_EPIMES/"
#reads_path="/media/marcgabriel/SAMSUNG/Dominika_XRN1_26_samples_2018/"
#reads_path="/media/marcgabriel/homeborn/Dominika_miRNA_seq/"
#reads_path="/media/marcgabriel/homeborn/Dominika_snoRNA_seq/"
#reads_path="/media/marcgabriel/saylar4/Nouritza_TNBClncRNA/"
#reads_path="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/"
#reads_path="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/"
#reads_path="/media/marcgabriel/saylar2/Nouritza_human_cell_lines/ /media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/"


#For FFPE files
#reads_path="/media/marcgabriel/saylar4/urines_ffpe_files/fastq_files/"

#For urine files
#reads_path="/media/marcgabriel/saylar4/urines_files/"

#reads_path="/media/marcgabriel/saylar5/prostate_cell_lines/"

#reads_path="/media/marcgabriel/saylar4/urines_ffpe_files/fastq_files/ /media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/"

#reads_path="/media/marcgabriel/saylar7/Hsieh_AC_GSE35469_PC3_cell_lines/trimmed_reads/"

#reads_path="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/urinary_EVs_clean_files_11_06_2021/urinary_EVs_fastq_files/"

#reads_path="/media/marcgabriel/saylar5/tissues_to_use_as_control_for_urines/"
#reads_path="/media/marcgabriel/saylar5/nuclear_cyto_prostate_samples_SRP199212/"
#reads_path="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/"
#reads_path="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/"
#reads_path="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/"
#reads_path="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021"
#reads_path="/media/marcgabriel/saylar8/plasticity_analysis/fastq_files/"
#reads_path="/media/marcgabriel/saylar8/plasticity_analysis/fastq_files_polyA/"
#reads_path="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/"
#reads_path="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/"
#reads_path="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/urinary_EVs_fastq_files/"
#reads_path="/media/marcgabriel/saylar8/Nouritza_EV_analysis/fastq_files/"
#reads_path="/media/marcgabriel/saylar9/Nouritza_breast_cancer_cell_lines_P_mestdagh_total_and_polyA/"
#reads_path="/media/marcgabriel/saylar9/Nouritza_human_TNBClncRNA_validation_batchPart2/" #these reads are in /media/marcgabriel/saylar10/Nouritza_TNBC_all_fastq_batches/ now

#reads_path="/media/marcgabriel/saylar11/Nouritza_TNBCdocPDX/fastq_files/"

#reads_path="/media/marcgabriel/saylar10/Nouritza_TNBC_cell_lines/"

#reads path for GTex alone
#reads_path="/media/marcgabriel/saylar15/GTEx_trimmed_fastq_611_samples_and_up_to_50/ /media/marcgabriel/saylar14/GTex_missing_DG_samples/449_missing_DG_samples/ "

#reads_path="/media/marcgabriel/saylar11/Rocco_RIPseq_analysis/fastq_files/"

#reads_path="/media/marcgabriel/saylar9/urines_files_new_batches_04_04_2023/ /media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b491/urinary_EVs_fastq_files/"

#reads path for GTex and Rocco lncap data
reads_path="/media/marcgabriel/saylar15/GTEx_trimmed_fastq_611_samples_and_up_to_50/ /media/marcgabriel/saylar14/GTex_missing_DG_samples/449_missing_DG_samples/ /media/marcgabriel/saylar5/LNCaP_data/"

#genome
#/home/marcgabriel/Documents/gencode19/gencode19_only_official_chromosomes.fa
#genome_fasta="/home/marcgabriel/Desktop/GFP.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode27lift37/all_genes.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode26/only_official_chr/genecode26_genes_fasta/all_genes.fa"
#genome_fasta="/media/marcgabriel/Transcend/Dominika_smallRNAseq_adapted_annotation/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/combined_annotations/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/combined_annotations/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/Nouritza_TNBC_human_pre_R_vs_pre_PD_dekupk_results/all_genes_exons_and_introns__contigs_pre_R_vs_pre_PD__pre_PD_vs_post_PD.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode32/gencode32_genes_fasta_format/all_genes.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode32/gencode32_genes_fasta_format/all_genes.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode32/gencode32_genes_fasta_format/gencode32_genes.fa"
#genome_fasta="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/gencode32_RepeatMasker.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode32/gencode.v32.transcripts_shortened_IDs.fa.gz"
#genome_fasta="/media/marcgabriel/saylar4/Nouritza_combined_annotation_gencode32_lncipedia_fantomcat_mitranscriptome_human_TNBC_scallop/combined_annotation_gencode32_lncipedia_mitranscriptome_fantom_scallop_spikesERCC_4456740.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode26/gencode.v26.transcripts_shortened_name.fa"
#genome_fasta="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/featurecounts_output_with_scallop_annot/DESeq2_results_with_scallop_annot_19_02_2021/intersection_of_comparisons/cyto/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/fasta_file_gencode26_DIS3_scallop/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar4/Nouritza_combined_annotation_gencode32_lncipedia_fantomcat_mitranscriptome_human_TNBC_scallop/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar7/Dominika_DIS3_analysis_22_03_2021/file.fa"
#genome_fasta="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/fasta_for_kallisto/all_genes.fa"
#genome_fasta="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/gencode26_DIS3_scallop_fasta_full_gene/all_genes.fa"

#genome_fasta="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/ORFs_analysis/TransDecoder_output/all_genes.facomplete_ORFs_transdecoder_genomic_coor_refined_CDS.fa"
#genome_fasta="/home/marcgabriel/Documents/gencode32/gencode32_genes_fasta_format/gencode32_genes.fa"

#genome_fasta="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/gencode26_DIS3_scallop_metatranscripts/gencode26_DIS3_scallop_metatranscripts_introns_numbered/intron_numberd_for_kallisto/all_genes.fa"

#genome_fasta="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt_extendedStartCodon/ribotricer_gencode26_DIS3_scallop_HCT116_eRNAs_conservative_refined_candidate_orfs.fa"

#genome_fasta="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/Dominika_riboseq_DIS3_diagenode_merged_runs_15_11_2021_and_16_12_2021/Dominika_ribotricer_analysis_gencode_scallop_eRNA_selected_length_from_unique_25_34nt_extendedStartCodon/ORFs_in_fasta_format/gencode26_scallop_eRNA_ORFs_extendedStartCodon.fa"


#genome_fasta="/media/marcgabriel/saylar8/Rocco_cut_and_run_analysis/gencode32_scallop_LNCaP_complete_27102021_transcripts.fa"

#genome_fasta="/media/marcgabriel/saylar5/LNCaP_data/gencode32_genes_and_PROCA11.fa"

genome_fasta="/media/marcgabriel/saylar5/LNCaP_data/fasta_from_LNCaP_complete_27102021/all_genes.fa"


#output files
#counter_outputs="/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/MAT1_overexp_EPIMES/GFP_analysis"
#counter_outputs="/media/marcgabriel/9d38aaa1-ce5d-4422-950f-0f083b797d13/oneTestAmongOthers"
#counter_outputs="/media/marcgabriel/homeborn/Dominika_miRNA_seq/kallisto_counts/"
#counter_outputs="/media/marcgabriel/homeborn/Dominika_snoRNA_seq/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar4/Nouritza_TNBClncRNA/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/kallisto_counts2/"
#counter_outputs="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/kallisto_counts2/"
#counter_outputs="/media/marcgabriel/saylar2/Nouritza_human_TNBClncRNA/kallisto_counts2/"
#counter_outputs="/media/marcgabriel/saylar2/Nouritza_human_cell_lines/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar2/Nouritza_human_genes_contigs_pre_R_vs_pre_PD__pre_PD_vs_post_PD/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode_and_repeats/"
#counter_outputs="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/"
#counter_outputs="/media/marcgabriel/saylar5/prostate_cell_lines/kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar5/tissues_to_use_as_control_for_urines/kallisto_counts/"

#counter_outputs="/media/marcgabriel/saylar5/tissues_to_use_as_control_for_urines/kallisto_counts_gencode32_RepeatMasker/"
#counter_outputs="/media/marcgabriel/saylar5/nuclear_cyto_prostate_samples_SRP199212/prostate_cyto_nuc_kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar5/nuclear_cyto_prostate_samples_SRP199212/prostate_cyto_nuc_kallisto_counts2/"
#counter_outputs="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/kallisto_counts_batch1And2_22_01_2021/"
#counter_outputs="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/featurecounts_output_with_scallop_annot/DESeq2_results_with_scallop_annot_19_02_2021/intersection_of_comparisons/cyto/scallop_final_common_transcripts_cyto_kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/fasta_file_gencode26_DIS3_scallop/"
#counter_outputs="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/gencode32_repeats_kallisto_counts/"
#counter_outputs="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/kallisto_counts_batch1And2_22_01_2021_2/"
#counter_outputs="/media/marcgabriel/saylar7/Dominika_DIS3_analysis_22_03_2021/gencode_and_scallop_transcripts_counts/"/media/marcgabriel/saylar9/Nouritza_breast_cancer_cell_lines_P_mestdagh_total_and_polyA/
#counter_outputs="/media/marcgabriel/saylar7/Anna_EVs_prostate_cell_lines_05_02_2021/kallisto_counts_gencode32_repeats_PNT2_PC3/"
#counter_outputs="/media/marcgabriel/saylar7/Nouritza_human_cell_lines/kallisto_counts_gencode_scallop/"
#counter_outputs="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/fasta_for_kallisto/"
#counter_outputs="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/fasta_for_kallisto_full_gene/"
#counter_outputs="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/kallisto_counts_gencode32_transcripts/ORFs_analysis/kallisto_counts_riboseq_ORFs/"
#counter_outputs="/media/marcgabriel/dd42f9e0-f4a1-4604-a0dc-c8768c723b49/urinary_EVs_clean_files_11_06_2021/kallisto_counts_6normaluEVs_6tumoruEVssamples/"
#counter_outputs="/media/marcgabriel/saylar5/Dominika_DIS3_analysis_10_12_2020/Dominika_DIS3_scallop_output3/cuffmerge_results/gencode26_DIS3_scallop_metatranscripts/gencode26_DIS3_scallop_metatranscripts_introns_numbered/intron_numberd_for_kallisto/kallisto_index_separate_introns_of_metatranscript/"
#counter_outputs="/media/marcgabriel/saylar8/plasticity_analysis/kallisto_counts2/"

#counter_outputs="/media/marcgabriel/saylar8/plasticity_analysis/kallisto_counts_early_rec_vs_prolonged_surv/"
#counter_outputs="/media/marcgabriel/saylar8/plasticity_analysis/kallisto_counts_polyA/"
#counter_outputs="/media/marcgabriel/saylar2/FFPE_vs_Urines_DESEq2/genes_expression2/kallisto_counts_on_transcripts"

#counter_outputs="/media/marcgabriel/saylar9/circRNAs_20_uEV_samples_JEV_revision/kallisto_counts_gencode_repeats/"

#counter_outputs="/media/marcgabriel/saylar8/Nouritza_EV_analysis/kallisto_counts_EVs/"

#counter_outputs="/media/marcgabriel/saylar9/Nouritza_breast_cancer_cell_lines_P_mestdagh_total_and_polyA/kallisto_counts/"

#counter_outputs="/media/marcgabriel/saylar9/Nouritza_human_TNBClncRNA_validation_batchPart2/kallisto_counts/"

#counter_outputs="/media/marcgabriel/saylar11/Nouritza_TNBCdocPDX/kallisto_counts/"

#counter_outputs="/media/marcgabriel/saylar10/Nouritza_TNBC_cell_lines/kallisto_counts/"

#counter_outputs="/media/marcgabriel/saylar17/Dominika_GTEX_kallisto_counts_ribotricer_gencode26_DIS3_scallop_HCT116_eRNAs_conservative_refined_candidate_orfs/"

#counter_outputs="/media/marcgabriel/saylar15/Dominika_GTEX_kallisto_counts_ribotricer_gencode26_DIS3_scallop_HCT116_eRNAs_conservative/"

#counter_outputs="/media/marcgabriel/saylar11/Rocco_RIPseq_analysis/Rocco_RIP_seq_kallisto_transcripts_counts/"

#counter_outputs="/media/marcgabriel/saylar9/Anna_uEVs_healthy_and_tumor_kallisto_counts/"

counter_outputs="/media/marcgabriel/saylar5/LNCaP_data/kallisto_GTex_and_LNCap_counts/"





chrom_size="/home/marcgabriel/Documents/gencode32/chrom_size.tsv"

#input_gff="/media/marcgabriel/saylar4/Nouritza_PDX_TNBClncRNA/combined_annotations/combined_annotation_gencode32_lncipedia_mitranscriptome_fantom.gff"
#input_gff="/media/marcgabriel/saylar4/Nouritza_human_TNBClncRNA/Nouritza_TNBC_human_pre_R_vs_pre_PD_dekupk_results/combined_annotation_Genelvl__contigs_pre_R_vs_pre_PD__pre_PD_vs_post_PD.gff"

input_gff="/home/marcgabriel/Documents/gencode32/gencode.v32.annotation_sorted.gff3"

#input_gff="/media/marcgabriel/saylar4/Nouritza_combined_annotation_gencode32_lncipedia_fantomcat_mitranscriptome_human_TNBC_scallop/combined_annotation_gencode32_lncipedia_mitranscriptome_fantom.gff"

#input_gff="/media/marcgabriel/saylar8/Rocco_cut_and_run_analysis/gencode32_scallop_LNCaP_complete_27102021.gff"


#e.g : pre_R_1_, post_PD_1_
#separator_after_file_name="_"
separator_after_file_name="_"

#separator_after_file_name="."

#number of process (samples) to run in parallel
process_number_limit=8

#threads used by kallisto for one sample
kallisto_threads=4

#kallisto program
#kallisto="/home/marcgabriel/new_dekupl/dekupl-run/bin/kallisto"
kallisto="/home/marcgabriel/Downloads/kallisto/kallisto"

gffread="/home/marcgabriel/Desktop/gffread-0.9.12.Linux_x86_64/gffread"

do_bam="no"

#kmer size for kallisto index (default=31)
#kmer_size=15
kmer_size=31

#library type ; if it's unstranded, just put "", otherwise choose between --rf-stranded & --fr-stranded
#--rf-stranded = mate2 is in the same direction as the RNA
#library_type="--fr-stranded"
library_type=""
#library_type="--rf-stranded"

#clean the existing bam file, if you think it's not in a standard way ("0" pos fort start or end, no string in the CIGAR == Ving gives eerros becaause of that)
clean_existing_bam="no"

check_file="no"
#check_file="yes"

#check_reads="yes"
check_reads="no"

#fragment length (for single-end)
fragment_length=100

sd_of_frag_len=80

###############  process the inputs ###################
#######################################################

#standardize the output dir name
counter_outputs="${counter_outputs}/"
counter_outputs=$(echo $counter_outputs |sed 's/\/\//\//g')
if [ ! -d $counter_outputs ]; then mkdir -p $counter_outputs; fi 

#directory of scripts for each sample (to be laucnhed in parallel)
sample_scripts_dir="${counter_outputs}/sample_scripts/"
if [ -d $sample_scripts_dir ]; then rm -rf $sample_scripts_dir; fi 
mkdir $sample_scripts_dir


######### dev ###########

#if no kallisto gtf, construct it
myGtf=${counter_outputs}kallisto_gtf.gtf

if [[ ! -f ${counter_outputs}kallisto_gtf.gtf ]];then

        #take only exons to facilitate the downstream processing
		grep -P "\texon\t" $input_gff >${counter_outputs}raw_gtf.gtf

		#gff to gtf
		$gffread ${counter_outputs}raw_gtf.gtf --force-exons -T -F -o ${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf.gtf

		cat ${counter_outputs}raw_gtf.gtf| sed 's/[;]/\t/g'|cut -f1-11 |sort -k1,1 -k10,10 -k11,11 -k7,7 | bedtools groupby -g 1,10,11,7 -c 4,5 -o min,max|awk -F"\t" 'OFS="\t"{print $1,".","gene",$5,$6,".",$4,".",$2";"" "$3";"}'|sed 's/ gene_id/gene_id/g' >${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf_genes.gtf 

		cat ${counter_outputs}raw_gtf.gtf| sed 's/[;]/\t/g' |cut -f1-11|sort -k1,1 -k10,10 -k11,11 -k9,9 -k7,7 | bedtools groupby -g 1,10,11,9,7 -c 4,5 -o min,max|awk -F"\t" 'OFS="\t"{print $1,".","transcript",$6,$7,".",$5,".",$2";"" "$3";"" "$4";"}' |sed 's/ gene_id/gene_id/g'>${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf_transcripts.gtf 

		cat ${counter_outputs}raw_gtf_genes.gtf ${counter_outputs}raw_gtf_transcripts.gtf ${counter_outputs}raw_gtf.gtf |sort -k1,1 -k4,4n -k4,5nr >${counter_outputs}raw_gtf.tmp && mv ${counter_outputs}raw_gtf.tmp ${counter_outputs}raw_gtf.gtf


		awk '( $3 ~ /gene/ )' ${counter_outputs}raw_gtf.gtf > $myGtf
		awk '( $3 ~ /transcript/ )' ${counter_outputs}raw_gtf.gtf >> $myGtf
		awk '( $3 ~ /exon/ && $7 ~ /+/ )' ${counter_outputs}raw_gtf.gtf | sort -k1,1 -k4,4n >> $myGtf
		awk '( $3 ~ /exon/ && $7 ~ /-/ )' ${counter_outputs}raw_gtf.gtf | sort -k1,1 -k4,4nr >> $myGtf


fi



##########"end of dev #######"

#look for all conds ; keep the order in the file
all_conds=($(grep -v "^#" $design|cut -f2 |awk '!seen[$0]++'))

echo -e "all conds are : ${all_conds[*]}\n"

#array to store all files
all_files=()

#array with all new names
all_new_names=()

#set -x
#for each sample of each cond, assign a name based on the condition (for conds normal & tumoral with 2 samples each : normal_1, normal_2...tumoral_1, tumoral_2)
for one_cond in $(seq 0 $((${#all_conds[*]}-1)));do

	#for a given cond, take all its files (samples)
	#files_one_cond=($(grep -P "\t${all_conds[$one_cond]}$" $design |cut -f1|sort -n))
	files_one_cond=($(grep -v "^#" $design |grep -P "\t${all_conds[$one_cond]}$" |cut -f1))
 
 
	#number of the 1st rep of the condition
	cond_start=1
 
	#number of the last rep of the condition
	cond_end=${#files_one_cond[*]}
 
	#store all the new names for this condition
	cond_new_names=()
 
    #loop across the samples of the cond, the for each of them concatenate the name of the condition & the number
	for ((i=$cond_start; i<=$cond_end; i++)); do cond_new_names+=(${all_conds[$one_cond]}_${i}) ; done
 
	all_files+=(${files_one_cond[*]})
 
	all_new_names+=(${cond_new_names[*]})

done

#echo -e "all files are : ${all_files[*]}"

#echo -e "new names are : ${all_new_names[*]}"

echo -e "new design : "

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g')

paste <(echo -e ${all_files[*]}|sed 's/ /\n/g') <(echo -e ${all_new_names[*]} |sed 's/ /\n/g') >${counter_outputs}new_design.tsv



echo -e "sample\tnb_reads_R1" >${counter_outputs}reads_summary.tsv



index_name=$(basename $genome_fasta|sed 's/\.fa//g'|sed 's/\.fasta//g')

#if the kallisto index of the genome isn't built, do it
index="${counter_outputs}${index_name}_kallisto_index.idx"

#build the index file if it's not there
if [ ! -f $index ];then 

   echo -e "\nwe are going to make the kallisto index with $genome_fasta\n..."

  $kallisto index -i $index --kmer-size=$kmer_size $genome_fasta

fi


#summary of the count files
>${counter_outputs}summary.tsv


##############  run a loop across the samples in order to process them #################
########################################################################################

all_results=()


#create a script for each sample, that will run Kallisto & process the output file
for one_sample in $(seq 0 $((${#all_files[*]}-1)));do

          #for a given sample of a given cond, take its new name
		  one_sample_name=$(echo "${all_new_names[$one_sample]}")
		  
		  #if [[ "$one_sample_name" != "Colon_49" ]];then
		  
		  #continue
		  
		  #fi
		  
		  library_type=""
		  
		  if [[ "$one_sample_name" =~ "LNCaP_FV" ]];then
		  
		    library_type="--rf-stranded"
		    
		    echo -e "we're going to change the lib type for LNcap Rocco...\n"
		  
		  fi
		  
		  #path to the sample
		  counter_outputs_OneSample="${counter_outputs}${one_sample_name}/"
		  
		  if [ ! -d $counter_outputs_OneSample ]; then mkdir $counter_outputs_OneSample; fi 
		  
		   echo -e "#!/bin/bash\n\n" >${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			#if the count file isn't there, use kallisto 
		   #echo -e "\nif [ ! -f ${counter_outputs_OneSample}abundance.tsv ] || [[ \$(wc -l ${counter_outputs_OneSample}abundance.tsv|awk '{print \$1}') -lt 1 ]];then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		   echo -e "\nif [ ! -f ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv ] || [[ \$(head ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv|wc -l|awk '{print \$1}') -lt 1 ]];then\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
		   
		   
		   
		   #set -x
		  		 
			  #fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}${separator_after_file_name}*" |grep ".*\.fastq"|sort -n)) || { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}${separator_after_file_name}*"|grep -v -i "report" |grep ".*\.fastq"|sort -n)) 
			  #|| { echo "no fastq files 1 !!" 1>&2; exit; }
			  
			  if [[ ${#fastq_list_OneSample[*]} -eq 0 ]];then
			  
			     echo -e "\n-> no fastq files with sep \"${separator_after_file_name}\" for $one_sample_name ($one_sample) , we're going to try with \".\"..."
			  
			     fastq_list_OneSample=($(find $reads_path -name "${all_files[$one_sample]}.*"|grep -v -i "report" |grep ".*\.fastq"|sort -n))|| { echo "no fastq files 2 !!" 1>&2; exit; }
			  
			  
			  fi
			  
			  #set +x
			  
			  echo -e "echo -e \"fastq files for condition $one_sample_name : \n ${fastq_list_OneSample[*]} \n--------------\n\"" >>${sample_scripts_dir}subscript_${one_sample_name}.sh

			  
			  #if [[ ${#fastq_list_OneSample[*]} -eq 2 ]];then
			  
			 
			  
				  #fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1\.fastq"|sort -u)
				  
				  
				  #set -x
				  #fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1_.*\.fastq"|sort -u|head -n1)
				  
				  fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "_R1.*\.fastq"|sort -u|head -n1)
				  
				  if [[ ! -f "$fastq_list_OneSample_R1" ]];then
				  
				    fastq_list_OneSample_R1=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "1\.fastq"|sort -u|head -n1)
				  
				  fi
				  
				   
				 
				  #fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2\.fastq"|sort -u)
				  
				  #fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2_.*\.fastq"|sort -u|head -n1)
				  
				  fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "_R2.*\.fastq"|sort -u|head -n1)
				  
				  
				  if [[ ! -f "$fastq_list_OneSample_R2" ]];then
				  
				    fastq_list_OneSample_R2=$(echo -e "${fastq_list_OneSample[*]}"|sed 's/ /\n/g'|grep "2\.fastq"|sort -u|head -n1)
				  
				  fi
				  #set +x
				  
				  if [[ "$check_file" == "yes" ]];then
				  
						if [[ $(gzip -t $fastq_list_OneSample_R1 2>&1|grep "not in") ]] || [[ $(gzip -t $fastq_list_OneSample_R1 2>&1|grep "unexpected end of") ]] ;then echo -e "-> issue with $fastq_list_OneSample_R1 !\n";fi
				  
				  fi
				  
				  if [[ "$check_reads" == "yes" ]];then
				  
						nb_reads=$(zcat $fastq_list_OneSample_R1|wc -l |awk '{OFMT="%f";print $1/4}' )
						
						echo -e "${one_sample_name}\t${nb_reads}" >>${counter_outputs}reads_summary.tsv
				  
				  fi
				  
	      
			      if [[ "$do_bam" == "yes" ]];then
			      
				  #run the quantif
			      echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --chromosomes=$chrom_size --genomebam --gtf=$myGtf -o ${counter_outputs_OneSample} ${fastq_list_OneSample_R1} ${fastq_list_OneSample_R2}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
			      #remove reads without start & end positions, or CIGAR string
			      echo -e "mv ${counter_outputs_OneSample}/pseudoalignments.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && rm ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh			      
			      
			      else
			      
			      
				  #run the quantif
			      echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads -o ${counter_outputs_OneSample} ${fastq_list_OneSample_R1} ${fastq_list_OneSample_R2}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			      
			      			      
			      
			      
			      
			      fi
				  
			  #elif [[ ${#fastq_list_OneSample[*]} -eq 1 ]];then
			  
			  
				  #if [[ "$check_file" == "yes" ]];then
				  
						#if [[ $(gzip -t ${fastq_list_OneSample[*]} 2>&1|grep "not in") ]];then echo -e "-> issue with $fastq_list_OneSample_R1 !\n";fi
				  
				  #fi			     
			  

			       
			       #if [[ "$do_bam" == "yes" ]];then
			       
			       #echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --chromosomes=$chrom_size --genomebam --gtf=$myGtf --single --fragment-length=$fragment_length --sd=$sd_of_frag_len -o ${counter_outputs_OneSample} ${fastq_list_OneSample[*]}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			       
			       #echo -e "mv ${counter_outputs_OneSample}/pseudoalignments.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && rm ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			       
			       
			       #else
			       
			       
			       #echo -e "${kallisto} quant -i $index $library_type -t $kallisto_threads --single --fragment-length=$fragment_length --sd=$sd_of_frag_len -o ${counter_outputs_OneSample} ${fastq_list_OneSample[*]}|| { echo \"kallisto failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh			       
			       
			       
			       #fi
			   
			  
			  #elif [[ ${#fastq_list_OneSample[*]} -eq 0 ]];then
			  
			     #echo -e "no fastq file for ${all_files[$one_sample]} $one_sample_name"
			     
			     ##exit 1
			     
			     #echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			     
			     ##it could be usefull if you just want to process some files in the list
			     #continue
			  
			  
			  #elif [[ ${#fastq_list_OneSample[*]} -gt 2 ]];then
			  
			    #echo -e "too many fastq files for ${all_files[$one_sample]} / $one_sample_name : ${fastq_list_OneSample[*]}"
			    
			    #exit 1
			  
			  #fi
			  
			  echo -e "cp ${counter_outputs_OneSample}abundance.tsv ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "\nfi\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh

			  
			  
			  
			  #process the count file in order to have columns like this : ID	counts
			  #echo -e "awk 'OFS=\"\\\t\"{if(NR>1){print \$1,\$4}}' ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv >${counter_outputs_OneSample}${one_sample_name}_abundance.tmp && mv ${counter_outputs_OneSample}${one_sample_name}_abundance.tmp ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  
			  echo -e "if [[ \$(head -n1 ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv|awk '{print NF}') -gt 2 ]];then\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  			  
				echo -e "awk 'OFS=\"\\\t\"{if(NR>1){print \$1,\$4}}' ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv >${counter_outputs_OneSample}${one_sample_name}_abundance.tmp && mv ${counter_outputs_OneSample}${one_sample_name}_abundance.tmp ${counter_outputs_OneSample}${one_sample_name}_abundance.tsv" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "fi\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  all_results+=(${counter_outputs_OneSample}${one_sample_name}_abundance.tsv)
			  
			  
				
			  #give rights to the script	
			  chmod 755 ${sample_scripts_dir}subscript_${one_sample_name}.sh
			  
			  echo -e "${counter_outputs_OneSample}${one_sample_name}_abundance.tsv\t${one_sample_name}" >>${counter_outputs}summary.tsv
			  
			  #ensure that the bam is well formed
			  if [[ "$do_bam" == "yes" ]];then
			  
			    #force to clean the existing bam file
			    if [[ "$clean_existing_bam" == "yes" ]];then
			  
			  			       echo -e "samtools view -@ $kallisto_threads -h ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam|awk 'OFS=\"\\\t\"{if((\$4!=0 && \$8 !=0 && \$6!=\"*\")||\$1~/^@/){print}}'|samtools view -@ $kallisto_threads -Sbh - >${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.test.bam && mv ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.test.bam ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam && samtools index ${counter_outputs_OneSample}/${one_sample_name}_pseudoalignments_sorted.bam || { echo \"kallisto bam failure for ${one_sample_name} !!\" 1>&2; exit; }\n\n" >>${sample_scripts_dir}subscript_${one_sample_name}.sh
			  			       
			    fi
			  
			  fi
					   
done



#run the subscripts in parallel

#here, if you are on a cluster, and you want to use a qsub command, just use a loop on the content of ${sample_scripts_dir}, and give the result to qsub
find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash sorting and/or counting subsritpts failure !" 1>&2; exit; }


all_IDs=($(cut -f2 ${counter_outputs}new_design.tsv))

echo -e "all samples are : ${all_IDs[*]}"





header=$(echo -e "${all_IDs[*]}"|sed 's/ /\t/g')

#final combined counts
echo -e "feature\t${header}" >${counter_outputs}header.tsv


#initialize the pos of the count files
#a=1

#header_table=(feature)

#for i in ${all_IDs[*]};do


    #one_file=$(find ${counter_outputs} -name "${i}_abundance.tsv"|grep -E "\/${i}_abundance.tsv")
    
    #echo -e "sample is : ${i}, file is $one_file"
    
    #if [[ -f $one_file ]];then


		#if [[ $a -eq 1 ]];then

			#cat $one_file >${counter_outputs}combined_counts.tsv
			
		#else
		
			#paste ${counter_outputs}combined_counts.tsv <(cut -f2 $one_file) >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv
		
		
		#fi
		
		#header_table+=(${i})
	  
	  
	  #a=$((a+1))
  
  
  #else
  
    #echo -e "${i}_abundance.tsv doesn't exist, next !\n"
  
  #fi

#done

#with too many files, you have to increase the limits
ulimit -Sn 10240 # The soft limit


#for a given line, from the 2nd column, if the column number has 0 after modulo 2, store the result (the count is at the columns with even numbers)
cat ${counter_outputs}header.tsv <(paste ${all_results[*]} | awk '{all=$1;for(i=2;i<=NF;i++){if(i%2==0){all=all"\t"$i}};print all}') >${counter_outputs}combined_counts.tsv



#cat ${counter_outputs}header.tsv ${counter_outputs}combined_counts.tsv >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv

#cat <(echo -e "${header_table[*]}"|sed 's/ /\t/g') ${counter_outputs}combined_counts.tsv >${counter_outputs}combined_counts.tmp && mv ${counter_outputs}combined_counts.tmp ${counter_outputs}combined_counts.tsv

exit


if [[ $(grep "\.gz" $genome_fasta|wc -l) -eq 0 ]];then


	cat <(echo -e "id\tlength") <(awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $1"\t"length($2)}}' $genome_fasta) >${counter_outputs}matching_id_length.tsv

else

	cat <(echo -e "id\tlength") <(zcat $genome_fasta|awk 'BEGIN{RS=">"}''{if(NR>1){sub("\n","\t"); gsub("\n",""); print $1"\t"length($2)}}') >${counter_outputs}matching_id_length.tsv
	
fi

echo -e "\n== check file : ${counter_outputs}combined_counts.tsv ==\n"






