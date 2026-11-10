# Week 11 — Multi-Omics Paper Analysis

**Zhou Hanyi · SUAT24000202**

The assignment allows any three of the four supplied papers. I answered Papers 1–3 and cite the relevant main figures, analytical methods, and results for each. The original Word submission is kept alongside this Markdown version: [`Week11_Homework_Short_Answers_with_Figures.docx`](Week11_Homework_Short_Answers_with_Figures.docx).

## Paper 1 — A spatiotemporal proteomic map of human adipogenesis

*Nature Metabolism* (2024), [DOI 10.1038/s42255-024-01025-8](https://doi.org/10.1038/s42255-024-01025-8)

### Q1. How were the data generated and processed?

Tryptic peptides were measured by label-free DIA mass spectrometry as 1-hour single shots across four adipocyte models (SGBS, hAPC, TERT-hAPC, hWA) and primary samples from seven donors, with biological triplicates. Spectronaut v15.7 in directDIA mode searched the UniProt human database at 1% q-value. This produced 5,979–7,061 protein groups per cell model, 3,638–6,403 in primary samples, and 8,268 in the merged dataset. QC looked acceptable: the median coefficient of variation was 12.3% and replicate Pearson correlation was 0.97. In Perseus, values were log2-transformed, missing values imputed, and significant proteins called by ANOVA at permutation FDR < 0.01 (Fig. 1a–e; Fig. 1a).

### Q2. What did the temporal analysis show?

PCA ordered samples along the differentiation trajectory (PC1 24.7%, PC2 13.8%), with lipid-droplet proteins driving PC1 (Fig. 1c,d). Of 6,569 significant proteins, 1,979 had conserved profiles—positive correlation in every pairwise model comparison. Hierarchical clustering split these into early, intermediate, and late groups. Fisher enrichment (BH FDR < 0.1) showed cell-cycle shutdown early, mitochondria and fatty-acid metabolism in the middle, and cholesterol plus lipid-droplet proteins late (Fig. 1g,h). The roughly 1,500 unchanged proteins were interpreted as relatively stable housekeeping machinery.

### Q3. What did transcriptome and spatial analysis add?

Protein–mRNA correlation was only 0.39–0.48 at individual time points, and 16% of profiles were negatively correlated, so mRNA alone cannot fully predict protein abundance. For localization, PCP fractionation plus SVM prediction mapped 5,530 proteins in adipocytes and 5,675 in preadipocytes (mean accuracy 96% in adipocytes and 88% in preadipocytes), and 1,323 proteins changed compartment during differentiation (Figs. 2–3).

For C19orf12, localization and interaction evidence placed it around lipid droplets and in contact with the TOM40/TOM22/TOM5 mitochondrial import machinery. Knockdown increased lipid-droplet number, reduced basal lipolysis by about half, and lowered maximal respiration. Etomoxir removed the respiratory difference, while palmitate and isoprenaline tests argued against a substrate-uptake or lipase defect. Together these experiments support a role in mitochondrial long-chain fatty-acid utilization (Fig. 6). In clinical cohorts, its expression also correlated inversely with body fat and HOMA-IR.

## Paper 2 — A universal gut-microbiome-derived signature predicts cirrhosis

*Cell Metabolism* (2020), [DOI 10.1016/j.cmet.2020.06.005](https://doi.org/10.1016/j.cmet.2020.06.005)

### Q1. How were the data generated and processed?

The UCSD cohort contained 163 participants, including 27 NAFLD-cirrhosis patients, 54 non-NAFLD controls, and first-degree relatives. Stool DNA was shotgun sequenced on MiSeq and processed through HUMAnN2: MetaPhlAn2 for taxonomy and DIAMOND/UniRef90 plus MetaCyc for function. This yielded 310 species and 7,355 functional pathways. Fecal metabolites were extracted, run by UPLC-qTOF MS, and processed in MS-DIAL with MS-Finder annotation, yielding 435 metabolites (Fig. 1).

### Q2. How was the predictive signature built and validated?

DESeq2 differential abundance plus random-forest feature selection (R `caret`, mean decrease in Gini, 10-fold cross-validation) selected 19 species and 17 metabolites; each signature reached training AUC 0.91. The 19 species plus age reached AUC 0.88 in relatives, 0.86 in the Chinese cohort, and 0.89 in the Italian cohort. Adding albumin raised the relative/Chinese models to 0.91/0.95. In the mixed-fibrosis cohort, the model separated cirrhosis from NAFL and fibrosis at 0.85/0.84, improving to 0.94/0.91 with AST (Fig. 5).

### Q3. What other analyses supported the finding?

HUMAnN2 identified altered microbial pathways, and RV-coefficient analysis tested concordance across cohorts (Fig. 2). Spearman correlation linked the 19 species and 17 metabolites in a network of 33 nodes and 338 edges (Fig. 4). A culture experiment added mechanistic support: increasing chenodeoxycholic acid caused *Ruminococcus gnavus* to produce more C18-sphingosine and 1H-indole-3-carboxaldehyde. Clinical variables (age, albumin, AST) improved the random-forest models. The associations remain predictive rather than causal.

## Paper 3 — Shared and disease-specific host gene–microbiome associations

*Nature Microbiology* (2022), [DOI 10.1038/s41564-022-01121-z](https://doi.org/10.1038/s41564-022-01121-z)

### Q1. How were the data generated and processed?

The study analyzed 208 paired mucosal samples: 88 CRC, 78 IBD, and 42 IBS. RNA-seq reads were trimmed, mapped with HISAT2 to GRCh38, counted with featureCounts, and DESeq2-VST transformed, retaining 12,513/11,985/12,429 genes in CRC/IBD/IBS. 16S OTUs were clustered in USEARCH/QIIME, contaminants removed, and counts CLR-transformed, leaving 235/121/238 taxa. Procrustes analysis found global host–microbe correspondence only in CRC (permutation P = 0.0001); Mantel testing was consistent (Fig. 1).

### Q2. How were associations identified?

Sparse CCA (R package PMA, LOOCV-tuned, BH FDR < 0.1) found group-level host–microbe co-variation components. Gene-wise LASSO (`glmnet`, LOOCV lambda), desparsified-lasso p-values, and stability selection produced 755, 1,295, and 441 associations at FDR < 0.1 in CRC, IBD, and IBS. Fisher enrichment identified disease-specific pathways, including NOTCH1 signaling in CRC, thrombin/PAR signaling in IBD, and arachidonic-acid metabolism in IBS (Figs. 2–3).

### Q3. What was shared and what was disease-specific?

Three shared host pathways (including oxidative phosphorylation and RAC1 signaling) and three shared taxa (*Peptostreptococcaceae*, *Streptococcus*, *Staphylococcus*) recurred across diseases. The same RAC1 pathway paired with oral taxa in CRC, *Granulicatella* in IBD, and *Bifidobacterium/Bacteroides* in IBS; other shared host genes clustered in VEGF, complement/coagulation, and cytokine pathways. Thus, a shared host immune-metabolic axis appears to be connected to disease-specific microbes. Because the study is cross-sectional and 16S does not resolve microbial function, directionality and functional mechanism remain unresolved.

## Appendix: cited figures

### Paper 1

![Paper 1 Fig. 1](figures/paper1_fig1_proteomic_dataset.png)

Fig. 1 | Proteomic dataset of four adipocyte models and primary samples.

![Paper 1 Fig. 2](figures/paper1_fig2_spatial_workflow.png)

Fig. 2 | Spatial proteomics workflow and organelle prediction (PCP–SVM).

![Paper 1 Fig. 3](figures/paper1_fig3_organelle_relocalization.png)

Fig. 3 | Organelle relocalization of 1,323 proteins during differentiation.

![Paper 1 Fig. 4](figures/paper1_fig4_bcaa_one_carbon.png)

Fig. 4 | Relocalization of BCAA and one-carbon metabolism enzymes.

![Paper 1 Fig. 5](figures/paper1_fig5_contact_sites_ld.png)

Fig. 5 | Mitochondrial contact-site proteins and lipid-droplet proteins.

![Paper 1 Fig. 6](figures/paper1_fig6_c19orf12.png)

Fig. 6 | C19orf12 characterization, knockdown phenotypes, and clinical cohorts.

### Paper 2

![Paper 2 Fig. 1](figures/paper2_fig1_cohort_overall.png)

Fig. 1 | Cohort design and overall microbiome–metabolome analysis.

![Paper 2 Fig. 2](figures/paper2_fig2_pathways_concordance.png)

Fig. 2 | Discriminatory microbial pathways and cross-cohort concordance.

![Paper 2 Fig. 3](figures/paper2_fig3_metabolites_diversity.png)

Fig. 3 | Discriminatory metabolites and diversity associations.

![Paper 2 Fig. 4](figures/paper2_fig4_network_rgnavus.png)

Fig. 4 | Correlation network and *R. gnavus* culture experiment.

![Paper 2 Fig. 5](figures/paper2_fig5_rf_validation.png)

Fig. 5 | Random-forest signature and external validation (AUC).

### Paper 3

![Paper 3 Fig. 1](figures/paper3_fig1_design_procrustes.png)

Fig. 1 | Study design and Procrustes analysis of paired samples.

![Paper 3 Fig. 2](figures/paper3_fig2_rac1_pathways.png)

Fig. 2 | Shared RAC1 pathway and disease-specific pathways.

![Paper 3 Fig. 3](figures/paper3_fig3_association_heatmaps.png)

Fig. 3 | Association heatmaps and disease-specific pathway networks.

![Paper 3 Fig. 4](figures/paper3_fig4_shared_microbes_genes.png)

Fig. 4 | Shared microbes and shared host genes across the three diseases.

## Submission record

- Repository: [santianwan/Bioinformatics_homework_zhouhanyi](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi)
- Commit: [`a11f7cc`](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi/commit/a11f7cc9870724a4baf7adca32563d4b78225246)
