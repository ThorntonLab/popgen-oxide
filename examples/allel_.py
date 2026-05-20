import allel

vcf = allel.read_vcf("ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
gt = allel.GenotypeArray(vcf['calldata/GT'])
ac = gt.count_alleles()
print(allel.mean_pairwise_difference(ac).sum())