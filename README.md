# vcfgl [![test](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml/badge.svg)](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml)


vcfgl is a lightweight utility tool for simulating genotype likelihoods from genotypes in VCF/BCF files. 

Given an input file with genotypes, average per-site read depth, and error rate parameters, the program can simulate the sampling of bases with base-calling errors and the genotype likelihoods.

You can print the genotype likelihoods in various formats:

- **GL tag**: Genotype likelihoods (enabled by default)
- **GP tag**: Genotype posterior probabilities (requires `-addGP 1`)
- **PL tag**: Phred-scaled genotype likelihoods (requires `-addPL 1`)

## Other features

- You can simulate unobserved invariable sites using `-explode 1`.

- You can define different depths for individuals using `--depths-file <filename>`. 
 
- Use `--depth inf` to set the values to known true variables.

- You can remove the unobserved alt alleles using `--trim-alt-alleles 1`.

- Simulate **QS tag**: Normalized phred-score allele quality sums (requires `-addQS 1`)

- Simulate **I16 tag**: Auxiliary tag used for calling (requires `-addI16 1`)


## Installation

```
git clone git@github.com:isinaltinkaya/vcfgl.git; cd vcfgl; make
```

## Usage


```
Usage:	vcfgl -i <input> [options]


Options:
	-v/--verbose		Verbose (default:0=disabled)
	-@/--threads		Number of threads (default:1)

    ->	Input/Output
	-i/--input		Input file (required)
	-o/--output		Output file prefix (default:output)
	-O/--output-mode	Output mode (default:b)
				v	VCF file
				b	BCF file
				z	Compressed VCF file (vcf.gz)
				u	Uncompressed BCF file

    ->	Simulation parameters
	-d/--depth		Mean per-site read depth
				Use `--depth inf` to set the simulated values to known true variables
	-df/--depths-file	File containing mean per-site read depth for each sample (conflicts with -d)
	-e/--error-rate		Error rate
	--error-qs		Simulate errors in quality scores (default:0=disabled)
				1	Simulate errors in quality scores by simulating the errors in the
					probability of wrong base calls not accounted in the reported quality 
					scores and genotype likelihoods. 
				2	Simulate errors in quality scores by simulating the errors in the 
					reported quality scores and genotype likelihoods. 
	--beta-variance		Variance of the beta distribution (default:disabled)
	--precise-gl		Should the program use the precise error rates in the  genotype likelihood
				 calculation? (default:1=enabled)
				1	Use precise error rates in the genotype likelihood calculation.
				0	Use the error rates calculated from the discretised phred scaled 
					quality scores in the genotype likelihood calculation.
	--pos0			Are the input coordinates are 0-based? (default:0=no)
				If input cordinates are 0 based, use --pos0 1 to shift positions by +1
	--seed			Random seed used to initialize the random number generator
	--platform		Simulate base qualities for a specific sequencing platform (default:0=disabled)
				1	NovaSeq 6000. The qualities are binned into four possible quality values: 2, 12, 23 and 37.
	-explode		Also simulate sites that were not observed in the input file (default:0=disabled)


Commands:
    ->	Specify which fields to simulate
	-addGL			Add GL field (default:1=enabled)
	-addGP			Add GP field (default:0=disabled)
	-addPL			Add PL field (default:0=disabled)
	-addI1			Add I16 field (default:0=disabled)
	-addQS			Add QS field (default:0=disabled)
	-addFormatDP		Add FORMAT/DP field (default:1=enabled)
	-addFormatAD		Add FORMAT/AD field (default:0=disabled)
	-addInfoDP		Add INFO/DP field (default:0=disabled)
	-addInfoAD		Add INFO/AD field (default:0=disabled)

```



## Input file

Input is typically a VCF/BCF file obtained by `tskit.write_vcf`.

Input file shoud have binary haplotypes set as REF and ALT alleles. For obtaining this format from mutation simulations, you can use the binary model with `sim_mutations`:

```python
mut_ts = msprime.sim_mutations(ts, ..., model="binary")
with open("simulated.vcf","w") as vcf:
    mut_ts .write_vcf(vcf)
```

For more information about the python code please refer to <https://tskit.dev/tskit/docs/stable/python-api.html>.


- **Input file position shifting**

Position shifting is controlled by `--pos0` option. Depending on your simulation setup, the VCF file generated tskit the positions may be 0-based (in contrast to regular vcf specs that is 1-based). If `-pos0 1` is used, all positions in the file will be shifted by `+1`.


## Simulate unobserved invariable sites (`-explode 1`)

You can use `-explode 1` to expand simulations to the sites not observed in the input VCF file.


Example command:
```
./vcfgl -i example.vcf -e 0.01 -d 5 -O v -d 3 --seed 42 -explode 1 
```

Example input: `example.vcf` 

Our input file contains two sites, and its header contains one chromosome named "chr22" with length 5.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=5>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
chr22	1	.	0	1	.	PASS	.	GT	0|0	0|0
chr22	3	.	0	1	.	PASS	.	GT	1|1	1|1
```


Example output: `output.vcf`

The missing 3 sites (at position 2,4, and 5) are generated by `-explode 1` by setting all individuals' genotypes to `0|0`. Thus the output file contains information for all sites.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=5>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##fileDate=Tue Sep 12 15:46:49 2023
##source=vcfgl --input example.vcf --output output --output-mode v --error-rate 0.010000 --depth 3.000000 --depths-file (null) --pos0 0 --seed 42 -explode 1 -printBaseCounts 0 -addGP 0 -addPL 0 -addI16 0 -addQS 0
##source=vcfgl version: v0.3.0 ee51b23-dirty
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
chr22	1	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551	0|0:5:0,-1.49785,-12.3638,-1.49785,-12.3638,-12.3638,-1.49785,-12.3638,-12.3638,-12.3638
chr22	2	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:3:0,-0.898711,-7.41827,-0.898711,-7.41827,-7.41827,-0.898711,-7.41827,-7.41827,-7.41827	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551
chr22	3	.	A	C,G,T	.	PASS	.	GT:DP:GL	1|1:1:-2.47276,-0.29957,0,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276	1|1:4:-9.89103,-1.19828,0,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103
chr22	4	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:5:0,-1.49785,-12.3638,-1.49785,-12.3638,-12.3638,-1.49785,-12.3638,-12.3638,-12.3638	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551
chr22	5	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:4:0,-1.19828,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103,-9.89103	0|0:4:0,-1.19828,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103,-9.89103
```


## Define different depths for individuals using `--depths-file <filename>`
 
Example: You want to set the average per-site read depth of different individuals to the target read depths as:

sample | target depth
-- | --
ind1 | 1 
ind2 | 5
ind3 | 0.4
ind4 | 4

Your depths file should contain one line per individual.

Example depths file: `depths.txt`
```
1
5
0.4
4
```

Then, you can run vcfgl using 
```
./vcfgl -i input.vcf --depths-file depths.txt
```


## Set the values to known true variables using `--depth inf`

This is useful for comparing the genotype calling-based methods with genotype likelihood-based methods. It will set the values in GL, GP, and PL tags to the best possible values corresponding to the input genotypes. Therefore, instead of using `--depth 100`, you can simulate the best possible values quickly.

___

## Genotype likelihood calculation

This program uses the same functions as in [ANGSD msToGlf](https://github.com/ANGSD/angsd/blob/master/misc/msToGlf.c) and 

The program is designed to produce simulated files equivalent to those simulated using [ANGSD msToGlf](https://github.com/ANGSD/angsd/blob/master/misc/msToGlf.c). Simulated genotype likelihoods follow the GATK definition.  

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

## Errors in quality scores

You can use `--error-qs 1` to simulate errors in quality scores by simulating the errors in the probability of wrong base calls not accounted in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution-sampled error rates and the reported quality scores use the beta distribution mean (i.e. the error rate).

You can use `--error-qs 2` to simulate errors in quality scores by simulating the errors in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution mean (i.e. the error rate) and the reported quality scores use the beta distribution-sampled error rates.


## Simulate base qualities for a specific sequencing platform

You can use `--platform 1` to simulate base qualities for a specific sequencing platform. Currently, only NovaSeq 6000 is supported. The qualities are binned into four possible quality values: 2, 12, 23 and 37.


## Example

```
./vcfgl -i test/t1.vcf -o test/t1_vcf_vcfgl -e 0.01 --seed 42 --depth 1 --pos0 1
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

MS equivalent:

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

### - Ancestral and derived alleles

By default, `vcfgl` treats the allelic state `0` in input file as `A`, and `1` as `C`. Therefore a genotype of `00` corresponds to `AA`, `01` to `AC`, `10` to `CA`, and `11` to `CC`.

Binary representation | Allele
-- | --
0 | A
1 | C
