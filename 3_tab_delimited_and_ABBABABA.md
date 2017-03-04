# Make a tab delimited file
```
~/tabix-0.2.6/bgzip XLXVXGXM_merged_sorted.bam.vcf.gz
~/tabix-0.2.6/tabix -f -p vcf XLXVXGXM_merged_sorted.bam.vcf.gz
zcat XLXVXGXM_merged_sorted.bam.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > XLXVXGXM_merged_sorted.bam.vcf.gz.tab
```
