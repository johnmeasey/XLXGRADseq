# here is some information on where the data are and where the XL genome is.

The raw data from the 2016 XLXG RAD plate is in this directory on info:
`/net/infofile4-inside/volume1/scratch/ben/2017_XL_XG_RADseq`

To demultiplex the data we used the `process_radtags` module of `Stacks`.

`/usr/local/bin/process_radtags -f 801_801_S316_L007_R1_001.fastq.gz -i gzfastq -b XLXG_barcodes -o ./samples/ -e sbfI -r -c -q --filter_illumina`

Using this barcode file (but without the name column`
```
TCCGGAGCGCTGCAGG	XL_CPT1
CTAACACGGCTGCAGG	XL_CPT1
AGCTTCGATTTGCAGG	XL_CPT1
TCGCCGCAATTGCAGG	XL_CTP2
TCAGTTCCGGTGCAGG	XL_CPT2
CGGAAGTGAGTGCAGG	XL_CPT2
GTTGCTAGACTGCAGG	XL_CPT3
AATAGATTCATGCAGG	XL_CPT3
AGCTGATACATGCAGG	XL_CPT3
ATCAGTAGAATGCAGG	XL_CPT4
TCGTCTTAGTTGCAGG	XL_CPT4
GCTCAGCCAGTGCAGG	XL_CPT4
CGGCTACTTCTGCAGG	XGUAE_59
CAAGCCGGTTTGCAGG	XGUAE_59
TTGCGCAAGCTGCAGG	XGUAE_59
TACGATGGAGTGCAGG	XGUAE_65
GCAATATACATGCAGG	XGUAE_65
AAGAATTCGGTGCAGG	XGUAE_65
TCGGCAGTCGTGCAGG	XGUAE_97
AGTTCCATTGTGCAGG	XGUAE_97
TTCTTGCGCTTGCAGG	XGUAE_97
AGCAATCTAATGCAGG	XGUAE_93
GAATTGTCGCTGCAGG	XGUAE_93
CTTCGACATATGCAGG	XGUAE_93
GAGATATGGTTGCAGG	XGUAE_92
CTCCTTGGAGTGCAGG	XGUAE_92
GTGTCTCTTGTGCAGG	XGUAE_92
TGCAGTTATCTGCAGG	XGUAE_72
TTCTGGAATATGCAGG	XGUAE_72
ACGCAACACATGCAGG	XGUAE_72
ACTGCCTCAATGCAGG	XGUAE_70
ACATCAATATTGCAGG	XGUAE_70
CCTCTTATCATGCAGG	XGUAE_70
TATCGTTAGTTGCAGG	XGUAE_71
TAGTGCGGTCTGCAGG	XGUAE_71
GGCCGGTAACTGCAGG	XGUAE_71
AGGAACCTCGTGCAGG	XGUAE_43
TTATCCGTAGTGCAGG	XGUAE_43
CGCTATACGGTGCAGG	XGUAE_43
CACGCAACGATGCAGG	XGUAE_42
TGTCCTAGGATGCAGG	XGUAE_42
ATCCGTCTACTGCAGG	XGUAE_42
GGACTCACGGTGCAGG	XGUAE_44
GCGTCCTGCCTGCAGG	XGUAE_44
ACTTGACCGGTGCAGG	XGUAE_44
AATGGTGACTTGCAGG	XGUAE_36
CTAACAGTATTGCAGG	XGUAE_36
TCATAGGCTATGCAGG	XGUAE_36
GCTGCACGGTTGCAGG	XG12_07
GCCGCAATGCTGCAGG	XG12_07
CGCTTCTCTGTGCAGG	XG12_07
CTCATTAACCTGCAGG	XG153
TTGATGGTGCTGCAGG	XG153
CAACATGAAGTGCAGG	XG153
ATGAAGGCAGTGCAGG	XG92
GATGGACTAATGCAGG	XG92
GCATGGAGGTTGCAGG	XG92
GTATATCCACTGCAGG	XGL713_177
CGTACCTTGCTGCAGG	XGL713_177
TATCGCGGAGTGCAGG	XGL713_177
CATGCATACTTGCAGG	XGL713_123
TCACTGAGAATGCAGG	XGL713_123
ATCCATAAGATGCAGG	XGL713_123
CTGTTAGATTTGCAGG	XGL713_123
TAACTGGTACTGCAGG	XGL713_179
GCCGGTGATTTGCAGG	XGL713_179
AGACGAATAGTGCAGG	XGL713_179
GGTCATTGTATGCAGG	XGL713_180
CGGTCGTTACTGCAGG	XGL713_180
GACGGACAGGTGCAGG	XGL713_180
ATTGCCACCGTGCAGG	XGL713_180
CAGAACCAGCTGCAGG	XGL713_181
GCTGTGCAGATGCAGG	XGL713_181
AAGACCAATCTGCAGG	XGL713_181
CGCGCGGCTGTGCAGG	XGUAE_124
GCATGAGGCGTGCAGG	XGUAE_124
TGACGGTGATTGCAGG	XGUAE_124
CTTACCGGAGTGCAGG	XLJONK_14
ACACACATCATGCAGG	XLJONK_14
TGTTGTCCGCTGCAGG	XLJONK_14
ATCCGCGACGTGCAGG	XLJONK_14
AAGTCGAGTATGCAGG	BJE3639
CCAATCAAGATGCAGG	BJE3639
TGAGCCAGCTTGCAGG	BJE3639
TCTAGAGAAGTGCAGG	BJE3545
GCCGAGGTGATGCAGG	BJE3545
GGTGAGTCGGTGCAGG	BJE3545
CTTATTCTACTGCAGG	XGL713_232
CAAGAGACGTTGCAGG	XGL713_232
GCACGTCTCCTGCAGG	XGL713_232
GTGCTCTCTATGCAGG	XGL713_232
GCGTCAGATGTGCAGG	XM_1
AATGAATCAGTGCAGG	XM_1
ATTAGGAGGCTGCAGG	XM_1
TTCTTCAGACTGCAGG	XM_1
CGCACTTGATTGCAGG	FGXCONTROL
```

# Preparing genome

Unzip the genome

`gunzip Xla.v91.repeatMasked.fa.gz`

and then index it 

`bwa index -a bwtsw Xla.v91.repeatMasked.fa`

`samtools faidx Xla.v91.repeatMasked.fa`

and now make a dict file:

`~/jre1.8.0_111/bin/java -jar /usr/local/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Xla.v91.repeatMasked.fa OUTPUT=Xla.v91.repeatMasked.dict`

# rename the demultiplexed files
```
mv sample_TCCGGAGCGC.fq XL_CPT1a.fq
mv sample_CTAACACGGC.fq XL_CPT1b.fq
mv sample_AGCTTCGATT.fq XL_CPT1c.fq
mv sample_TCGCCGCAAT.fq XL_CTP2a.fq
mv sample_TCAGTTCCGG.fq XL_CPT2b.fq
mv sample_CGGAAGTGAG.fq XL_CPT2c.fq
mv sample_GTTGCTAGAC.fq XL_CPT3a.fq
mv sample_AATAGATTCA.fq XL_CPT3b.fq
mv sample_AGCTGATACA.fq XL_CPT3c.fq
mv sample_ATCAGTAGAA.fq XL_CPT4a.fq
mv sample_TCGTCTTAGT.fq XL_CPT4b.fq
mv sample_GCTCAGCCAG.fq XL_CPT4c.fq
mv sample_CGGCTACTTC.fq XGUAE_59a.fq
mv sample_CAAGCCGGTT.fq XGUAE_59b.fq
mv sample_TTGCGCAAGC.fq XGUAE_59c.fq
mv sample_TACGATGGAG.fq XGUAE_65a.fq
mv sample_GCAATATACA.fq XGUAE_65b.fq
mv sample_AAGAATTCGG.fq XGUAE_65c.fq
mv sample_TCGGCAGTCG.fq XGUAE_97a.fq
mv sample_AGTTCCATTG.fq XGUAE_97b.fq
mv sample_TTCTTGCGCT.fq XGUAE_97c.fq
mv sample_AGCAATCTAA.fq XGUAE_93a.fq
mv sample_GAATTGTCGC.fq XGUAE_93b.fq
mv sample_CTTCGACATA.fq XGUAE_93c.fq
mv sample_GAGATATGGT.fq XGUAE_92a.fq
mv sample_CTCCTTGGAG.fq XGUAE_92b.fq
mv sample_GTGTCTCTTG.fq XGUAE_92c.fq
mv sample_TGCAGTTATC.fq XGUAE_72a.fq
mv sample_TTCTGGAATA.fq XGUAE_72b.fq
mv sample_ACGCAACACA.fq XGUAE_72c.fq
mv sample_ACTGCCTCAA.fq XGUAE_70a.fq
mv sample_ACATCAATAT.fq XGUAE_70b.fq
mv sample_CCTCTTATCA.fq XGUAE_70c.fq
mv sample_TATCGTTAGT.fq XGUAE_71a.fq
mv sample_TAGTGCGGTC.fq XGUAE_71b.fq
mv sample_GGCCGGTAAC.fq XGUAE_71c.fq
mv sample_AGGAACCTCG.fq XGUAE_43a.fq
mv sample_TTATCCGTAG.fq XGUAE_43b.fq
mv sample_CGCTATACGG.fq XGUAE_43c.fq
mv sample_CACGCAACGA.fq XGUAE_42a.fq
mv sample_TGTCCTAGGA.fq XGUAE_42b.fq
mv sample_ATCCGTCTAC.fq XGUAE_42c.fq
mv sample_GGACTCACGG.fq XGUAE_44a.fq
mv sample_GCGTCCTGCC.fq XGUAE_44b.fq
mv sample_ACTTGACCGG.fq XGUAE_44c.fq
mv sample_AATGGTGACT.fq XGUAE_36a.fq
mv sample_CTAACAGTAT.fq XGUAE_36b.fq
mv sample_TCATAGGCTA.fq XGUAE_36c.fq
mv sample_GCTGCACGGT.fq XG12_07a.fq
mv sample_GCCGCAATGC.fq XG12_07b.fq
mv sample_CGCTTCTCTG.fq XG12_07c.fq
mv sample_CTCATTAACC.fq XG153a.fq
mv sample_TTGATGGTGC.fq XG153b.fq
mv sample_CAACATGAAG.fq XG153c.fq
mv sample_ATGAAGGCAG.fq XG92a.fq
mv sample_GATGGACTAA.fq XG92b.fq
mv sample_GCATGGAGGT.fq XG92c.fq
mv sample_GTATATCCAC.fq XGL713_177a.fq
mv sample_CGTACCTTGC.fq XGL713_177b.fq
mv sample_TATCGCGGAG.fq XGL713_177c.fq
mv sample_CATGCATACT.fq XGL713_123a.fq
mv sample_TCACTGAGAA.fq XGL713_123b.fq
mv sample_ATCCATAAGA.fq XGL713_123c.fq
mv sample_CTGTTAGATT.fq XGL713_123d.fq
mv sample_TAACTGGTAC.fq XGL713_179a.fq
mv sample_GCCGGTGATT.fq XGL713_179b.fq
mv sample_AGACGAATAG.fq XGL713_179c.fq
mv sample_GGTCATTGTA.fq XGL713_180a.fq
mv sample_CGGTCGTTAC.fq XGL713_180b.fq
mv sample_GACGGACAGG.fq XGL713_180c.fq
mv sample_ATTGCCACCG.fq XGL713_180d.fq
mv sample_CAGAACCAGC.fq XGL713_181a.fq
mv sample_GCTGTGCAGA.fq XGL713_181b.fq
mv sample_AAGACCAATC.fq XGL713_181c.fq
mv sample_CGCGCGGCTG.fq XGUAE_124a.fq
mv sample_GCATGAGGCG.fq XGUAE_124b.fq
mv sample_TGACGGTGAT.fq XGUAE_124c.fq
mv sample_CTTACCGGAG.fq XLJONK_14a.fq
mv sample_ACACACATCA.fq XLJONK_14b.fq
mv sample_TGTTGTCCGC.fq XLJONK_14c.fq
mv sample_ATCCGCGACG.fq XLJONK_14d.fq
mv sample_AAGTCGAGTA.fq BJE3639a.fq
mv sample_CCAATCAAGA.fq BJE3639b.fq
mv sample_TGAGCCAGCT.fq BJE3639c.fq
mv sample_TCTAGAGAAG.fq BJE3545a.fq
mv sample_GCCGAGGTGA.fq BJE3545b.fq
mv sample_GGTGAGTCGG.fq BJE3545c.fq
mv sample_CTTATTCTAC.fq XGL713_232a.fq
mv sample_CAAGAGACGT.fq XGL713_232b.fq
mv sample_GCACGTCTCC.fq XGL713_232c.fq
mv sample_GTGCTCTCTA.fq XGL713_232d.fq
mv sample_GCGTCAGATG.fq XM_1a.fq
mv sample_AATGAATCAG.fq XM_1b.fq
mv sample_ATTAGGAGGC.fq XM_1c.fq
mv sample_TTCTTCAGAC.fq XM_1d.fq

```
concatenate the samples

```
cat XL_CPT1a.fq XL_CPT1b.fq XL_CPT1c.fq > concat/XL_CPT1.fq 
cat XL_CTP2a.fq XL_CPT2b.fq XL_CPT2c.fq > concat/XL_CPT2.fq 
cat XL_CPT3a.fq XL_CPT3b.fq XL_CPT3c.fq > concat/XL_CPT3.fq
cat XL_CPT4a.fq XL_CPT4b.fq XL_CPT4c.fq > concat/XL_CPT4.fq
cat XGUAE_59a.fq XGUAE_59b.fq XGUAE_59c.fq > concat/XGUAE_59.fq
cat XGUAE_65a.fq XGUAE_65b.fq XGUAE_65c.fq > concat/XGUAE_65.fq
cat XGUAE_97a.fq XGUAE_97b.fq XGUAE_97c.fq > concat/XGUAE_97.fq
cat XGUAE_93a.fq XGUAE_93b.fq XGUAE_93c.fq > concat/XGUAE_93.fq
cat XGUAE_92a.fq XGUAE_92b.fq XGUAE_92c.fq > concat/XGUAE_92.fq  
cat XGUAE_72a.fq XGUAE_72b.fq XGUAE_72c.fq > concat/XGUAE_72.fq
cat XGUAE_70a.fq XGUAE_70b.fq XGUAE_70c.fq > concat/XGUAE_70.fq
cat XGUAE_71a.fq XGUAE_71b.fq XGUAE_71c.fq > concat/XGUAE_71.fq
cat XGUAE_43a.fq XGUAE_43b.fq XGUAE_43c.fq > concat/XGUAE_43.fq
cat XGUAE_42a.fq XGUAE_42b.fq XGUAE_42c.fq > concat/XGUAE_42.fq
cat XGUAE_44a.fq XGUAE_44b.fq XGUAE_44c.fq > concat/XGUAE_44.fq
cat XGUAE_36a.fq XGUAE_36b.fq XGUAE_36c.fq > concat/XGUAE_36.fq
cat XG12_07a.fq XG12_07b.fq XG12_07c.fq > concat/XG12_07.fq
cat XG153a.fq XG153b.fq XG153c.fq > concat/XG153.fq
cat XG92a.fq XG92b.fq XG92c.fq > concat/XG92.fq
cat XGL713_177a.fq XGL713_177b.fq XGL713_177c.fq > concat/XGL713_177.fq
cat XGL713_123a.fq XGL713_123b.fq XGL713_123c.fq XGL713_123d.fq > concat/XGL713_123.fq
cat XGL713_179a.fq XGL713_179b.fq XGL713_179c.fq > concat/XGL713_179.fq
cat XGL713_180a.fq XGL713_180b.fq XGL713_180c.fq XGL713_180d.fq > concat/XGL713_180.fq
cat XGL713_181a.fq XGL713_181b.fq XGL713_181c.fq > concat/XGL713_181.fq
cat XGUAE_124a.fq XGUAE_124b.fq XGUAE_124c.fq > concat/XGUAE_124.fq
cat XLJONK_14a.fq XLJONK_14b.fq XLJONK_14c.fq XLJONK_14d.fq > concat/XLJONK_14.fq
cat BJE3639a.fq BJE3639b.fq BJE3639c.fq > concat/BJE3639.fq
cat BJE3545a.fq BJE3545b.fq BJE3545c.fq > concat/BJE3545.fq
cat XGL713_232a.fq XGL713_232b.fq XGL713_232c.fq XGL713_232d.fq > concat/XGL713_232.fq
cat XM_1a.fq XM_1b.fq XM_1c.fq XM_1d.fq > concat/XM_1.fq
```

# Align the fq files using bwa mem and sort them too

```
#!/bin/bash                                                                                            

path_to_data="samples/concat"
path_to_chromosome="/net/infofile4-inside/volume1/scratch/ben/2017_XL_XG_RADseq/XL_v9.1/"
chromosome="Xla.v91.repeatMasked"

individuals="XL_CPT1
XL_CPT2
XL_CPT3
XL_CPT4
XGUAE_59
XGUAE_65
XGUAE_97
XGUAE_93
XGUAE_92
XGUAE_72
XGUAE_70
XGUAE_71
XGUAE_43
XGUAE_42
XGUAE_44
XGUAE_36
XG12_07
XG153
XG92
XGL713_177
XGL713_123
XGL713_179
XGL713_180
XGL713_181
XGUAE_124
XLJONK_14
BJE3639
BJE3545
XGL713_232
XM_1
BJE1488
BJE1489
BJE261
BJE263
BJE264
BJE265
BJE266
BJE267"

for each_individual in $individuals
do

echo ${each_individual}
    bwa mem -M -t 16 -r "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}\tPL:illumina" $path_to_chromosome/$chromosome.fa $path_to_data/${each_individual}.fq | samtools view -bSh - > $path_to_data/${each_individual}.bam
    samtools sort $path_to_data/${each_individual}.bam $path_to_data/${each_individual}_sorted
    samtools index $path_to_data/${each_individual}_sorted.bam
done

```

# Make a merged bam file

`samtools merge XLXVXGXM_merged_sorted.bam *_sorted.bam` (see problems next!)

This did not really work for some reason.  The merged file had only one sample in it.  So we tried picard instead:

`~/jre1.8.0_111/bin/java -jar /usr/local/picard-tools-1.131/picard.jar MergeSamFiles I=BJE1488_sorted.bam I=BJE1489_sorted.bam I=BJE261_sorted.bam I=BJE263_sorted.bam I=BJE264_sorted.bam I=BJE265_sorted.bam I=BJE266_sorted.bam I=BJE267_sorted.bam I=BJE3545_sorted.bam I=BJE3639_sorted.bam I=XG12_07_sorted.bam I=XG153_sorted.bam I=XG92_sorted.bam I=XGL713_123_sorted.bam I=XGL713_177_sorted.bam I=XGL713_179_sorted.bam I=XGL713_180_sorted.bam I=XGL713_181_sorted.bam I=XGL713_232_sorted.bam I=XGUAE_124_sorted.bam I=XGUAE_36_sorted.bam I=XGUAE_42_sorted.bam I=XGUAE_43_sorted.bam I=XGUAE_44_sorted.bam I=XGUAE_59_sorted.bam I=XGUAE_65_sorted.bam I=XGUAE_70_sorted.bam I=XGUAE_71_sorted.bam I=XGUAE_72_sorted.bam I=XGUAE_92_sorted.bam I=XGUAE_93_sorted.bam I=XGUAE_97_sorted.bam I=XL_CPT1_sorted.bam I=XL_CPT2_sorted.bam I=XL_CPT3_sorted.bam I=XL_CPT4_sorted.bam I=XLJONK_14_sorted.bam I=XM_1_sorted.bam O=XLXVXGXM_merged_sorted.bam `


# index the bam file

`samtools index XLXVXGXM_merged_sorted.bam`

`~/samtools_2016/bin/samtools mpileup -d8000 -ugf /net/infofile4-inside/volume1/scratch/ben/2017_XL_XG_RADseq/XL_v9.1/Xla.v91.repeatMasked.fa -t DP,AD XLXVXGXM_merged_sorted.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o XLXVXGXM_merged_sorted.bam.vcf.gz`
