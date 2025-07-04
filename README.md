# South African Blood Regulatory (SABR) Resource

The South African Blood Regulatory (SABR) Resource generated whole genome sequencing (WGS) and blood RNA-seq data from over 600 individuals spanning three South Eastern Bantu-speaking groups. This was a collaboration between [Variant Bio](https://www.variantbio.com/) and Michele Ramsay's group at [Wits University](https://www.wits.ac.za/), and includes individuals from the [AWI-Gen](https://h3africa.org/index.php/awi-gen/) cohort. These data were used to map genetic variants that impact gene expression, splicing, and cell type levels. A full description of the resource can be found in our [manuscript](https://www.nature.com/articles/s41588-025-02223-0).

**A map of blood regulatory variation in South Africans enables GWAS interpretation. Castel _et al._ 2025. _Nature Genetics_.**

> Functional genomics resources are critical for interpreting human genetic studies, but currently they are predominantly from European-ancestry individuals. Here we present the South African Blood Regulatory (SABR) resource, a map of blood regulatory variation that includes three South Eastern Bantu-speaking groups. Using paired whole genome and blood transcriptome data from over 600 individuals, we map the genetic architecture of 40 blood cell traits derived from deconvolution analysis, as well as expression, splice, and cell type interaction quantitative trait loci. We comprehensively compare SABR to the Genotype Expression (GTEx) Project and characterize thousands of regulatory variants only observed in African-ancestry individuals. Finally, we demonstrate the increased utility of SABR for interpreting African-ancestry association studies by identifying putatively causal genes and molecular mechanisms through colocalization analysis of blood-relevant traits from the Pan-UK Biobank. Importantly, we make full SABR summary statistics publicly available to support the African genomics community.

## Table of Contents

1. [Cohort and Data Overview](#cohort-and-data-overview)
2. [Summary Statistics](#summary-statistics)
    - [Genotype Data](#genotype-data)
    - [GWAS](#gwas)
    - [QTL Mapping](#qtl-mapping)
3. [Analysis Results](#analysis-results)
    - [xCell Disease Modeling](#xcell-disease-modeling)
    - [Colocalization Analyses](#colocalization-analyses)
4. [Zenodo](#zenodo)
5. [Controlled Access Data](#controlled-access-data)
6. [Methods](#methods)
    - [Genotype Calling and Imputation](#genotype-calling-and-imputation)
    - [Expression Quantification, Splice Quantification, and QTL Mapping](#expression-quantification-splice-quantification-and-qtl-mapping)
    - [xCell GWAS](#xcell-gwas)
    - [Colocalization](#colocalization)

## Cohort and Data Overview

The SABR cohort consists of 754 individuals who participanted in the AWI-Gen cohort and were recontacted and reconsented for this study. Venous whole blood samples were taken and used to carry out WGS at a median depth of 5.1x and paired-end, stranded RNA-sequencing with globin and rRNA depletion at a median depth of 30M mapped read pairs. All data are made available for non-commercial use only.

**All coordinates are provided in GRCh38**

- Individual-level sequencing data quality control metrics and inclusion in downstream analyses - [GitHub](data_tables/S1_participant_metadata.txt)

## Summary Statistics

Variant-level summary statics from WGS genotyping and full summary statistics from xCell GWAS and QTL mapping are publicly available via the links listed below. All full summary statistics files (`.txt.bgz`) are provided with a corresponding index file (`.tbi`). Note, genome-wide summary files are provided through an AWS S3 bucket and will require the AWS CLI [tool](https://aws.amazon.com/cli/) to download.

### Genotype Data
- South African enriched, putatively functional alleles - [GitHub](data_tables/S2_functional_alleles.txt)
- Variant-level summary statistics from imputed mid-pass WGS including allele frequencies and functional annotations - `s3://public.us-prod.variantbio.com/SABR/VARS/SABR_variant_summary.txt.bgz`
- Description of variant-level summary statistics fields - [GitHub](data_tables/variant_summary_stats.txt)

### GWAS
- List of xCell types included in analyses - [GitHub](data_tables/S3_xcell_cell_types.txt)
- xCell codes - [GitHub](data_tables/xcell_codes.txt)
- Summary of genome-wide significant loci (p < 5e-8) identified - [GitHub](data_tables/S5_xcell_gwas_loci.txt)
- Full summary statistics outputted by Hail for each GWAS, including p-value, beta (alt allele), standard error, minor allele frequency - `s3://public.us-prod.variantbio.com/SABR/XCELL_GWAS/`

### QTL Mapping
Gene-level results and full variant level summary statistics outputted by [fastQTL](https://github.com/francois-a/fastqtl) are provided for all QTL mapping runs.

1. Expression QTLs (eQTLs)
    - Gene-level results - [GitHub](data_tables/S7_cis_eqtl_genes.txt)
    - eVariant annotations - [GitHub](data_tables/S8_cis_eqtl_variants.txt)
    - Conditionally independent eQTLs - [GitHub](data_tables/S12_ind_eqtls.txt)
    - Conditionally independent eVariant annotations - [GitHub](data_tables/S13_ind_eqtl_variants.txt)
    - Nominally significant structural variant eQTLs - [GitHub](data_tables/S11_cis_eqtl_sv_results.txt)
    - Full summary statistics - `s3://public.us-prod.variantbio.com/SABR/EQTL/SABR_eQTL_allpairs.txt.bgz`
    - Conditionally independent eQTL summary statistics - `s3://public.us-prod.variantbio.com/SABR/EQTL/SABR_eQTL_conditional_variants.txt.gz`
2. Splice QTLs (sQTLs)
    - Gene-level results - [GitHub](data_tables/S9_cis_sqtl_genes.txt)
    - sVariant annotations - [GitHub](data_tables/S10_cis_sqtl_variants.txt)
    - Full summary statistics - `s3://public.us-prod.variantbio.com/SABR/SQTL/SABR_sQTL_allpairs.txt.bgz`
3. Cell-type Interaction eQTLs (ieQTLs)
    - Gene-level results - [GitHub](data_tables/S14_cis_ieqtl_genes.txt)
    - ieVariant annotations - [GitHub](data_tables/S15_cis_ieqtl_variants.txt)
    - Full summary statistics - `s3://public.us-prod.variantbio.com/SABR/IEQTL/`

## Analysis Results

### xCell Disease Modeling
- Modeling results for each xCell type by disease - [GitHub](data_tables/S4_xcell_disease_modeling.txt)

### Colocalization Analyses

1. SABR QTLs x PAN-UKBB African GWAS
    - List of African GWAS included in analysis - [GitHub](data_tables/S16_panukbb_afr_gwas.txt)
    - Colocalization results - [GitHub](data_tables/S17_afr_coloc_results.txt)
    - Colocalization lead variant annotations - [GitHub](data_tables/S18_afr_coloc_variants.txt)
2. SABR eQTLs x PAN-UKBB Multi-ancestry GWAS
    - List of multi-ancestrty (MA) GWAS included in analysis - [GitHub](data_tables/S19_panukbb_ma_gwas.txt)
    - Colocalization results - [GitHub](data_tables/S20_ma_coloc_results.txt)

## Zenodo

In addition to being avilable here and on AWS, eQTL and sQTL summary statistics have been uploaded to a Zenodo [repository](https://zenodo.org/records/15334125).

## Controlled Access Data

Indvidual-level data are availabe to authorized users purusing reasearch in line with informed consent and ethical approvals. Data are provided through the European Phenome Genome Archive (EGA) project [EGAS50000001008](https://ega-archive.org/studies/EGAS50000001008). The following is a list of data availble via controlled access.

- Joint called and imputed genotype VCF
    - HAIL-ZAAG-WGS-MP-1-2S.metadata.txt
    - HAIL-ZAAG-WGS-MP-1-2S.vcf.bgz
    - HAIL-ZAAG-WGS-MP-1-2S.vcf.bgz.tbi
- Extended indivividual-level metadata, including disease status
    - SABR_subject_phenotypes.txt
- xCell enrichment scores
    - SABR_xcell_quantifications.txt
- Expression quantifications (counts, TMM, TPM)
    - SABR_gene_reads.gct.gz
    - SABR_gene_tmm.bed.gz
    - SABR_gene_tpm.gct.gz
- Splice quantifications (junction counts, clusters)
    - SABR_leafcutter_junctions.bed.gz
    - SABR_leafcutter_clusters.bed.gz
- xCell GWAS input files (covariates)
    - SABR_gwas_covariates.txt
- QTL mapping input files (covariates, normalized quantifications)
    - SABR_EQTL_combined_covariates.txt
    - SABR_EQTL_expression.bed.gz
    - SABR_SQTL_combined_covariates.txt
    - SABR_SQTL_leafcutter.bed.gz
- Whole-genome sequencing data (FASTQs)
- RNA-sequencing data (FASTQS)

## Methods

All methods are described in detail in the Supplementary Materials provided with our [manuscript](xxx). Below, we briefly describe the methods used and link to relevant software and pipeline pages.

### Genotype Calling and Imputation

Genotype calling and imputation from mid-pass WGS data was carried out as described in [Emde et al.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07949-9). Code for mid-pass genotype calling and imputation is available on [GitHub](https://github.com/variant-bio/mid-pass).

### Expression Quantification, Splice Quantification, and QTL Mapping

The GTEx/TOPMed v8 [pipeline](https://github.com/broadinstitute/gtex-pipeline) was used for read mapping, quantifying expression levels, normalizing data, and mapping eQTLs with [fastQTL](https://github.com/francois-a/fastqtl). RNA-SeQC (v2.3.6) was used for gene quantification of collapsed genes with GENCODE v34 annotations. Splicing quantification and sQTL mapping was carried out using the approach described by the [GTEx Consortium](https://www.science.org/doi/10.1126/science.aaz1776). Cell type estimation and interaction expression QTL mapping was carried using the approach described by [Kim-Hellmuth et al.](https://www.science.org/doi/10.1126/science.aaz8528). For eQTL and ieQTL mapping, 60 PEER factors and 20 genotype PCs were used as covariates in addition to age, sex, and mean WGS depth. For sQTL mapping, 15 PEER factors and 20 genotype PCs were used as covariates in addition to age, sex, and mean WGS depth.

### xCell GWAS

All cell types with enrichment scores > 0 in > 50% of participants were used for GWAS. Cell type enrichment scores were inverse normal transformed and GWAS were run using the `linear_regression_rows()` function in [Hail](https://hail.is/) v0.2 with the following covariates: age, sex, sex\*age, sex\*age^2, mean WGS depth, and 20 genotype PCs.

### Colocalization

A subset of African and multi-ancestry GWAS from the [Pan-UKBB](https://pan.ukbb.broadinstitute.org/) were used for colocalization analysis. Colocalization analysis was carried out using [coloc](https://chr1swallace.github.io/coloc/) v5 for all significant QTLs (FDR < 5%) within 500kb of a significant (p<5e-8, multi-ancesrtry GWAS) or suggestive (p<5e-6, African ancestry GWAS) association signal. Minor allele frequencies, p-values, and a prior12 value of 1e-5 were used when running coloc.