# Prado-et-al-Postglacial-colonization-by-Peromyscus-leucopus

Scripts and input files used in the paper submitted to the Journal of Mammalogy:

Prado et al. Postglacial colonization in the Great Lakes Region by the white-footed mouse (Peromyscus leucopus): conflicts between genomic and field data 

(1) Bioinformatics

1-CFF_Create_whitelist_for_stacks2e.R: Create a whitelist to remove SNP’s positioned at the end of all loci due to an artificially increased number of SNPs observed at these last positions. Also, it removes loci with high theta values (above the 95 percentiles), given these are suggestive of sequencing and assembly errors.

(2) Diversity&Struture

2-genetic_diversity.R: Check the genetic diversity summary statistics for significant differences, and plots their means values and standard errors.

3-PCAs_projected_axes.R: Perform Principal Component Analysis (PCA) projecting the individuals from the two expanded populations onto the axes of the historic populations, following Lipson et al. 2018 recomendations

4- Fst&Mantel.R: Perform Mantel test with a sequential population dropout procedure

(3) Expansion

5- Expansion.R: Calculate the directionality index ψ (Peter and Slatkin 2013). This statistic detects the allele frequency clines created by successive founder events, where the further away a population is from the origin of the range expansion, the higher the probability that a SNP increases in allele frequency or becomes fixed. 

(4) Divergence_Time

6-moving&best_lhood.R: Move files among folders to facilitate Fastsimcoal analysis

Input files used in the Fastsimcoal runs

