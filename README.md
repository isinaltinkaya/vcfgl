# vcfgl [![test](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml/badge.svg)](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml)
## Small program to simulate genotype likelihoods from VCF GT tags

## Installation

```
git clone git@github.com:isinaltinkaya/vcfgl.git; cd vcfgl; make
```

## Usage

```
./vcfgl -in {INPUT} -mode {OUTPUT_MODE} -out {OUTPUT_FILE_PREFIX} -err {ERROR_RATE} -depth {DEPTH} -seed {SEED} -isSim {0|1} 
```

Arguments:

* **INPUT**
 
Input is a VCF/BCF file obtained by `tskit.write_vcf`.

After simulating a Tree Sequence, save your tree sequence as a VCF file:
```
with open("simulated.vcf", "w") as vcf_file:
    tree_sequence.write_vcf(vcf_file, ploidy=2)
```

For more information please refer to https://tskit.dev/tskit/docs/stable/python-api.html.

This program uses the GT tags in `simulated.vcf` to simulate genotype likelihoods.


* **OUTPUT_FILE_PREFIX**

File prefix for the output and arguments file the program will generate. 

e.g. Using `-out myout` will generate the files `myout.vcf` and `myout.arg`.
e.g. Using `-out myout` with `-mode b` will generate the files `myout.bcf` and `myout.arg`.


* **OUTPUT_MODE**

 Output modes can be set by `-mode {b,u,v,z}`. By default, a compressed BCF file will be generated.

```
b       compressed BCF [default]
u       uncompressed BCF
v       uncompressed  VCF
z       compressed VCF (bgzf compressed)
```

* **ERROR_RATE**

Error rate to be used in simulation. For more details, see [Genotype likelihood calculation](#genotype-likelihood-calculation).

* **DEPTH**

Mean per site read depth to be used in simulations.

* **SEED**

Random seed used to initialize the random number generator.

* **isSim**

In the VCF file generated by tskit, the positions are 0-based, in contrast to regular vcf specs that is 1-based. If `isSim 1` is used, all positions in the file will be shifted by `+1`.


## Genotype likelihood calculation

This program uses the same functions as in [ANGSD msToGlf](https://github.com/ANGSD/angsd/blob/master/misc/msToGlf.c) and is designed to produce equivalent results.

 
Genotype likelihoods follow the GATK definition.  

![image](https://user-images.githubusercontent.com/33371500/170043998-dcff8c7d-b483-42e4-b312-38b280970fc8.png)

![image](https://user-images.githubusercontent.com/33371500/170044034-89650362-5e0d-4321-810b-698a7fba0d20.png)


The genotype likelihood simulation procedure is as follows:

1. For each site; for each haplotype; sample number of reads from Poisson distribution with mean equal to given depth value

2. Simulating error

(2.1) Probability of error: e

(2.2) Probability of sampling an incorrect allele due to error: error rate × 1/3

(2.3) Probability of sampling correct allele: 1 − error rate


4. Normalize the log10 genotype likelihood by substracting the maximum genotype likelihood value
observed at site.

## Example


```
./vcfgl -in test/t1.vcf -out test/t1_vcf_vcfgl -err 0.01 -seed 42 -depth 1 -isSim 1
```

Input file: 

```
$ cat test/t1.vcf
##fileformat=VCFv4.2
##source=tskit 0.4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=51229805>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ind
chr22   282     .       0       1       .       PASS    .       GT      0|0
chr22   310     .       0       1       .       PASS    .       GT      1|0
chr22   587     .       0       1       .       PASS    .       GT      0|1
chr22   773     .       0       1       .       PASS    .       GT      1|1
```

Output file: 

```
$ cat test/t1_vcf_vcfgl.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.4.1
##contig=<ID=chr22,length=51229805>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated read depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ind
chr22   283     .       A       C,G,T   .       PASS    .       GT:DP:GL        0|0:1:0,-0.29957,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276,-2.47276
chr22   311     .       A       C,G,T   .       PASS    .       GT:DP:GL        1|0:0:.,.,.,.,.,.,.,.,.,.
chr22   588     .       A       C,G,T   .       PASS    .       GT:DP:GL        0|1:2:-1.87362,0,-1.87362,-2.17319,-2.17319,-4.34637,-2.17319,-2.17319,-4.34637,-4.34637
chr22   774     .       A       C,G,T   .       PASS    .       GT:DP:GL        1|1:1:-2.47276,-0.29957,0,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276
```


MS file equivalent of the test case:
```
$ cat test/t1.ms
ms 2 1
0

//
segsites: 4
positions: 0.0000 0.0000 0.0000 0.0000
0101
0011
```

### Msprime - Mutation simulations

`vcfgl` currently only works with simulated vcf files that have binary haplotypes set as REF and ALT. For obtaining this format from mutation simulations, you can use the binary model with `sim_mutations`:

```
mutts = msprime.sim_mutations(ts, rate=msp_u, model="binary")
with open(vcf_file,"w") as vcf:
                                mutts.write_vcf(vcf)
```

### - Order of genotype likelihoods

The order of genotype likelihoods in VCF file format is calculated as in the following pseudocode:

```
for P=ploidy and N=number of alternate alleles;
for a_p in 0:N; for a_p-1 in 0:a_p; print(a1,a2);
```

For P=2 N=3

```
0,1,2,3
A,C,G,T
00,01,11,02,12,22,03,13,23,33
AA,AC,CC,AG,CG,GG,AT,CT,GT,TT
```
