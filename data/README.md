# _Xantusia_ data

- `data` > `input_files`
  - `mtDNA_v7.csv` : file with mtDNA samples matched to metadata (e.g., latitude, longitude, original paper)
  - `mtDNA.Xantusia_outgroups.aln.fasta`: file used to infer mtDNA gene tree
  - `snapp.Xantusia_speciestree.snps`: file used to infer SNAPP species tree
  - `xantusia_samples_v10.csv` file with ddRAD samples matched to metadata (e.g., latitude, longitude, OTU)
  - `Xantusia.miss0.5.no_outgroups.vcf.zip`: variants for ingroup samples only - this variant set was used for most analyses (PCA, admixture, dstat)
  - `Xantusia.miss0.5.with_outgroups.phy.zip`: concatenated locus file for ingroup & outgroup files - this was used for phylogenetic inference
- `data` > `output_files`
  - `alphahull.ranges`: range files for alphahulls of species
  - `mtDNA.Xantusia_outgroups.aln.treefile`: IQtree inferred gene tree for mtDNA alignment
  - `nDNA.Xantusia_outgroups.tre`: IQtree inferred nDNA tree for ddRAD concatenated alignment (includes outgroups)
  - `snapp.mcc.Xantusia_speciestree.tre`: SNAPP species tree, results combined across multiple SNAPP runs