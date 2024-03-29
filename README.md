# vcf2smc_pro.py

An improvement to Terhorst et al.'s [vcf2smc](https://github.com/popgenmethods/smcpp). It is done by assigning a more accurate ancestral allele when creating the input files for smc++. It also has an experimental setting to correct for low coverage of highly homozygous individuals by utilizing pseudo-haploidism.

To run this program, you will need 5 arguments, with an optional 6th. NOTE: This can only be done on one contig at a time. All files can only contain information about that contig.
* **Names List** - A .txt file documenting the distinguished individuals and populations
* **VCF File** - This is a file of all called SNPS on the individuals you want to run.
* **Ancestral BED** - A file of the VCF positions in the first three columns, and the LiftOver positions on the fourth column
* **Chromosome Converter** - A space-delimited file with corresponding chromosome names in the ancestral FASTA and VCF
* **Ancestral Fasta** - An ancestral sequence of the population in question
* **Mask** (optional) - A standard mask file in BED format

## Creating Input Data:

### Names List:
To create the names list file, put the two distinguished individuals (AS NAMED IN THE VCF) on the first line, space delimited, followed by their ploidy. In general, the ploidy should be set to d (Diploid) for now. The program also allows for using h (pseudohaploidy), but it is not recommended. On the following lines, make the same pattern, as shown below:
```
Dist_ind_1,d Dist_ind_2,d
Pop1_ind1,d Pop1_ind2,d
Pop2_ind1,d Pop2_ind2,d
```
### VCF File:
Input a VCF File with indels filtered out. Software that is capable of doing this is [samtools mpileup](http://www.htslib.org/doc/samtools.html),[GATK UnifiedGenotyper](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php), and [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

### Ancestral BED:
This file is a combination of two different files - a BED file of the positions of the VCF, and the liftOver positions on what would be the ancestral genome species. Here is a sample of how to make it:
```
1. [zcat | ]awk '{if($0!~"#") print $1"\t"$2-1"\t"$2"\t"$1"_"$2}' example.vcf[.gz] > SNPs.bed
```
NOTE: If your ancestral genome is the same species, you can stop here. The next liftOver step is only if the ancestral sequence is in the form of another species.
```
2. liftOver SNPs.bed map_to_anc_seq_species.chain Anc_bed.bed unmapped
```
Your final file should look something like this:
```
chr5    175532976       175532977       chr2_1170
chr5    175532921       175532922       chr2_1225
chr5    175532805       175532806       chr2_1341
chr5    175532786       175532787       chr2_1360
chr5    180111015       180111016       chr2_2259
chr5    180111260       180111261       chr2_2486
chr5    180111302       180111303       chr2_2528
chr5    180111472       180111473       chr2_2674
chr5    180111632       180111633       chr2_2828
chr5    180111826       180111827       chr2_3026
chr5    180111829       180111830       chr2_3029
chr5    180111871       180111872       chr2_3068
chr5    180111896       180111897       chr2_3093
chr5    180112043       180112044       chr2_3241
chr5    180112107       180112108       chr2_3305
chr5    180112233       180112234       chr2_3421
chr5    180112552       180112553       chr2_3746
chr5    180112589       180112590       chr2_3783
```

Of course, if your ancestral sequence is the same species, the chromosome in the first column should be the same as the chromosome in the 5th column, as should the number in the third column be equal to the number in the 5th column.

### Chromosome Converter:
A space delimited file where the first column are the chromosome names present in the VCF file, and the second column are the corresponding chromosome names in the ancestral FASTA file. In cases where the names are the same, there would be two identical columns:

```
chr1a chr1a
chr2 chr2
chr3 chr3
chr4 chr4
chr5 chr5
chr6 chr6
chr7b chr7b
chr8 chr8
chr9 chr9
chr10 chr10
chr11 chr11
chr12 chr12
chr13 chr13
chr14 chr14
chr15 chr15
chr16 chr16
chr17 chr17
chr18 chr18
chr19 chr19
chr20 chr20
chr21 chr21
chr22a chr22a
chr23 chr23
```

### Ancestral Fasta:
The inferred ancestral sequence. This can be done by utilizing pairwise alignments and creating a sort of "consensus sequence" between closely related species. In our example, we used the [hg19 multiz100way alignment](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/), taking only homologous regions between nomLeu3, hg19 and rheMac3. Sites where hg19 and nomLeu3 disagreed, rheMac3 was used to reorient them.

After the SMC files are created, you can move onto the next step, using [smc++](https://github.com/popgenmethods/smcpp).
