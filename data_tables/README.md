## Supplementary Data Tables

### Table S1. SABR Participant Metadata and Sequencing Stats
Metadata for SABR participants including age, sex, SEB group (group). WGS statistics including mean depth (wgs_mean_dp), proportion of genotype calls imputed (wgs_impute_rate), number of non-imputed genotype calls (wgs_non_imputed). RNA QC metrics including RIN (rna_RIN), percentage of fragments over 200bp (rna_DV200), and number of RNA-seq mapped read pairs (rna_mapped_read_pairs). For each participant, the analyses they were included in is defined (in_wgs = genotype only analyses, in_xcell = cell type deconvolution GWAS, in_qtl = in QTL mapping).

[S1_participant_metadata.txt](S1_participant_metadata.txt)

### Table S2. South African Functional Alleles
Annotated variants that are predicted to have high functional impact (CADD > 30), present at MAF > 1% in SABR and unobserved in 1000 Genomes high coverage whole genomes. Variants were queried in ClinVar by rsID for any associated entries, and any conditions and classifications are listed (clinvar_condition, clinvar_classification). Allele frequencies are provided for the alternative (ALT) allele.

[S2_functional_alleles.txt](S2_functional_alleles.txt)

### Table S3. xCell Cell Types
Cell type enrichments that are estimated by xCell, with types included in downstream SABR analyses indicated (included_in_analyses).

[S3_xcell_cell_types.txt](S3_xcell_cell_types.txt)

### Table S4. xCell Disease Modeling
Results from a logistic regression model of disease ~ cell type for each cell type by disease comparison. P-values were adjusted for multiple testing across all comparisons using the FDR method.

[S4_xcell_disease_modeling.txt](S4_xcell_disease_modeling.txt)

### Table S5. xCell GWAS Loci
All loci that were genome-wide significant in at least one xCell GWAS (p < 5e-8) with top GWAS variant by p-value (lead_variant_id) and corresponding annotations listed. For each association Open Targets was queried to determine if the locus contained any associations with blood cell counts (ot_locus_associations) and if the lead variant was associated with blood cell counts (ot_variant_associations). Allele frequencies and betas are provided for the alternative (ALT) allele.

[S5_xcell_gwas_loci.txt](S5_xcell_gwas_loci.txt)

### Table S6. cis-QTL Mapping Summary
For each QTL mapping run, summary of number of significant tested genes (tested_genes), QTL genes (FDR < 5%, sig_genes), unique lead QTL variants (unique_lead_variants), tested_genes with minor allele count >= 10 (tested_genes_MAC10), significant QTL genes with MAC >= 10 (sig_genes_MAC10), unique lead QTL variants with MAC >= 10 (unique_lead_vars_MAC10), and minor allele frequency (MAF) in each SEB group, or the entire SABR cohort (ALL).

[S6_cis_qtl_mapping_summary.txt](S6_cis_qtl_mapping_summary.txt)

### Table S7. cis-eQTL Results
Gene-level permutation results from eQTL mapping using fastQTL. Columns are the standard outputs from fastQTL run in permutation mode. In column names “ma” means “minor allele”. Slope is always provided for the alternative allele.

[S7_cis_eqtl_genes.txt](S7_cis_eqtl_genes.txt)

### Table S8. eVariant Annotations
Annotations for all unique, lead variants for significant eQTLs (eVariants). Variants present in eQTLGen full summary statistics by rsID lookup are indicated (eqtlgen). Allele frequencies are provided for the alternative (ALT) allele.

[S8_cis_eqtl_variants.txt](S8_cis_eqtl_variants.txt)

### Table S9. cis-sQTL Results
Gene-level permutation results from sQTL mapping using fastQTL. Columns are the standard outputs from fastQTL run in permutation mode. In column names “ma” means “minor allele”. Slope is always provided for the alternative allele.

[S9_cis_sqtl_genes.txt](S9_cis_sqtl_genes.txt)

### Table S10. sVariant Annotations
Annotations for all unique, lead variants for significant sQTLs (sVariants). Variants present in eQTLGen full summary statistics by rsID lookup are indicated (eqtlgen). Allele frequencies are provided for the alternative (ALT) allele.

[S10_cis_sqtl_variants.txt](S10_cis_sqtl_variants.txt)

### Table S11. cis-eQTL Structural Variant Results
Variant-level results from cis-eQTL mapping for nominally significant structural variants (SVs). Only eGenes (FDR < 5%) were included and only SVs with nominal p-value < nominal p-value threshold for the eGene as determined by permutation mapping are reported. Allele frequencies are provided for the alternative (INS/DEL) allele. The method used for genotype calling is indicated by the “sv_method” column: GATK or PanGenie (PANG).

[S11_cis_eqtl_sv_results.txt](S11_cis_eqtl_sv_results.txt)

### Table S12. Independent cis-eQTLs
Gene-level results from conditional eQTL analysis, with nominal p-value threshold used for conditional analysis (p_value_threshold), number of conditionally independent eQTLs mapped (independent_signals), and variant IDs (variant_ids). Includes a flag for if an error occurred causing conditional analysis to be stopped before hitting 20 independent eQTLs (error).

[S12_ind_eqtls.txt](S12_ind_eqtls.txt)

### Table S13. Independent eVariant Annotations
Annotations for all unique, lead variants for conditionally independent eQTLs. Variants present in eQTLGen full summary statistics by rsID lookup are indicated (eqtlgen). Allele frequencies are provided for the alternative (ALT) allele.

[S13_ind_eqtl_variants.txt](S13_ind_eqtl_variants.txt)

### Table S14. cis-ieQTL Results
Gene-level results from ieQTL mapping using fastQTL with multiple testing correction using eigenMT with number of independent tests estimated by eigenMT (num_eigenmt_tests) and ieQTL p-value corrected for this number of tests (pval_eigenmt). Other columns are the standard outputs from fastQTL run in nominal pass mode. In column names “ma” means “minor allele”. Slope is always provided for the alternative allele.

[S14_cis_ieqtl_genes.txt](S14_cis_ieqtl_genes.txt)

### Table S15. ieVariant Annotations
Annotations for all unique, lead variants for significant ieQTLs (ieVariants). Variants present in eQTLGen full summary statistics by rsID lookup are indicated (eqtlgen). Allele frequencies are provided for the alternative (ALT) allele.

[S15_cis_ieqtl_variants.txt](S15_cis_ieqtl_variants.txt)

### Table S16. PAN-UKBB African-Ancestry GWAS
Summary information for African GWAS used in colocalization analysis including trait category (phenotype_category), regression type (model), number of genome-wide significant hits (hits), number of loci with at least one genome-wide significant hit (loci), and full path to summary statistics (path_sumstats).

[S16_panukbb_afr_gwas.txt](S16_panukbb_afr_gwas.txt)

### Table S17. PAN-UKBB African-Ancestry Colocalization Results
Results with PP4 > 0.50 from colocalization analysis including variant explaining maximum proportion of PP4 (coloc_variant), proportion of H4 explained by coloc_variant (coloc_variant_H4), QTL beta (beta_qtl), GWAS beta (beta_gwas), QTL nominal p-value (pvalue_qtl), and GWAS p-value (pvalue_gwas). Betas are provided for the alternative allele (ALT).

[S17_afr_coloc_results.txt](S17_afr_coloc_results.txt)

### Table S18. PAN-UKBB African-Ancestry Colocalization Variant Annotations
Annotations for all unique coloc_variants from colocalizations with PP4 > 0.50. Allele frequencies are provided for the alternative (ALT) allele.

[S18_afr_coloc_variants.txt](S18_afr_coloc_variants.txt)

### Table S19. PAN-UKBB Multi-Ancestry GWAS
Summary information for multi-ancestry GWAS used in colocalization analysis including trait category (phenotype_category), regression type (model), number of genome-wide significant hits (hits), number of loci with at least one genome-wide significant hit (loci), and full path to summary statistics (path_sumstats).

[S19_panukbb_ma_gwas.txt](S19_panukbb_ma_gwas.txt)

### Table S20. PAN-UKBB Multi-Ancestry Colocalization Results
Results with PP4 > 0.50 from colocalization analysis including variant explaining maximum proportion of PP4 (coloc_variant), proportion of H4 explained by coloc_variant (coloc_variant_H4), QTL beta (beta_qtl), GWAS beta (beta_gwas), QTL nominal p-value (pvalue_qtl), and GWAS p-value (pvalue_gwas). Betas are provided for the alternative allele (ALT).

[S20_ma_coloc_results.txt](S20_ma_coloc_results.txt)