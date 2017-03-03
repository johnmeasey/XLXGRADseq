# Make a tab delimited file
```
~/tabix-0.2.6/bgzip nonrecal_filtered_chr11_final.vcf
~/tabix-0.2.6/tabix -f -p vcf temp.bam.vcf.gz
zcat temp.bam.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > temp.bam.vcf.gz.tab
```
