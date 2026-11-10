# Week 11 Homework — Papers 1–3

**Zhou Hanyi · SUAT24000202**

## Paper 1: A spatiotemporal proteomic map of human adipogenesis

**Q1. How were the data generated and processed?**

They used label-free DIA mass spectrometry on four adipocyte cell models and primary samples from seven donors. Spectronaut/directDIA searched the human UniProt database, and the final merged dataset had 8,268 protein groups. The QC looked reasonable: median CV about 12.3%, and replicate correlations around 0.97. Downstream analysis in Perseus used log2 values, imputation, and ANOVA with FDR < 0.01 (Fig. 1a–e).

**Q2. What did the temporal analysis show?**

PCA followed the differentiation process quite clearly: PC1 was 24.7% and separated early from late samples (Fig. 1c,d). Of 6,569 significant proteins, 1,979 had similar temporal profiles across models. Clustering split them into early, intermediate, and late groups. The early group was mostly cell-cycle proteins going down, the middle group was mitochondrial/fatty-acid metabolism, and the late group had cholesterol synthesis and lipid-droplet proteins (Fig. 1g,h).

**Q3. What did the transcriptome and spatial analysis add?**

The protein–mRNA correlation was only 0.39–0.48, and 16% of profiles were actually negative, so transcript levels alone do not explain protein abundance. Spatial proteomics by PCP/SVM mapped 5,530 proteins in adipocytes and 5,675 in preadipocytes, with mean accuracy 96%/88%; 1,323 proteins changed compartment during differentiation (Figs. 2–3). C19orf12 was the interesting one: it localized near lipid droplets and mitochondria, its knockdown changed lipid accumulation and basal lipolysis, and etomoxir removed the respiration difference. That points to a role in long-chain fatty-acid use rather than a general lipase defect (Fig. 6).

![Paper 1, Fig. 1](figures/paper1_fig1_proteomic_dataset.png)

## Paper 2: A universal gut-microbiome-derived signature predicts cirrhosis

**Q1. How were the data generated and processed?**

The UCSD cohort had 163 participants, including 27 NAFLD-cirrhosis patients, 54 controls, and relatives. Shotgun stool DNA was processed with HUMAnN2/MetaPhlAn2, giving 310 species and 7,355 pathways. Fecal metabolites were measured by UPLC-qTOF and annotated with MS-DIAL/MS-Finder, giving 435 metabolites (Fig. 1).

**Q2. How was the predictive model built and validated?**

They used DESeq2 for differential features and random forests in R's `caret` package, with 10-fold CV and Gini importance. The final model had 19 species and 17 metabolites, each reaching training AUC 0.91. The 19-species model plus age worked in external cohorts too: AUC 0.88 in relatives, 0.86 in the Chinese cohort, and 0.89 in the Italian cohort. Adding albumin or AST improved the discrimination further (Fig. 5).

**Q3. What else supported the result?**

The pathway results were concordant across cohorts by RV-coefficient analysis, and the species/metabolite signature formed a 33-node, 338-edge Spearman correlation network (Figs. 2 and 4). The culture experiment was the most convincing part for me: increasing chenodeoxycholic acid made *R. gnavus* produce more C18-sphingosine and 1H-indole-3-carboxaldehyde. The model is still predictive, not causal.

![Paper 2, Fig. 5](figures/paper2_fig5_rf_validation.png)

## Paper 3: Shared and disease-specific host gene–microbiome associations

**Q1. How were the data generated and processed?**

They collected 208 paired mucosal samples: 88 CRC, 78 IBD, and 42 IBS. RNA-seq was mapped with HISAT2, counted with featureCounts, and transformed with DESeq2-VST; 16S OTUs were clustered in USEARCH/QIIME and CLR-transformed. Procrustes analysis found a global host–microbe match only in CRC (P = 0.0001), not IBD or IBS (Fig. 1).

**Q2. How were the associations identified?**

Two models were used. Sparse CCA looked for groups of host genes and microbes that varied together. Gene-wise LASSO plus stability selection then gave 755, 1,295, and 441 associations in CRC, IBD, and IBS at FDR < 0.1. Enrichment analysis found 18 disease-specific pathways, including NOTCH1 signaling in CRC, thrombin/PAR signaling in IBD, and arachidonic-acid metabolism in IBS (Figs. 2–3).

**Q3. What was shared and what was disease-specific?**

Only a small part was shared. Three host pathways, including oxidative phosphorylation and RAC1 signaling, appeared in all three diseases, and three taxa—*Peptostreptococcaceae*, *Streptococcus*, and *Staphylococcus*—were also shared. The surprise was that the same RAC1 pathway connected to different microbes in each disease: oral taxa in CRC, *Granulicatella* in IBD, and *Bifidobacterium/Bacteroides* in IBS. So the host pathway can be shared while the microbial partners change. Since this is cross-sectional 16S data, it cannot prove direction or microbial function.

![Paper 3, Fig. 2](figures/paper3_fig2_rac1_pathways.png)

## Submission record

- Repository: [santianwan/Bioinformatics_homework_zhouhanyi](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi)
- Commit: [`a35fccd`](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi/commit/a35fccd305fb9254bde92f7a95737753ca89f36b)
