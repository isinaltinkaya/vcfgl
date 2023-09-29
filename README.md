# vcfgl [![test](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml/badge.svg)](https://github.com/isinaltinkaya/vcf-gl/actions/workflows/test.yml)


vcfgl is a lightweight utility tool for simulating genotype likelihoods from genotypes in VCF/BCF files. Given read depth (`--depth`) and error rate (`--error-rate`) parameters, the program can be used to simulate and print genotype likelihoods in various formats:

- genotype likelihoods (GL tag) (default)
- genotype posterior probabilities (GP tag) (requires `-addGP`)
- phred-scaled genotype likelihoods (PL tag) (requires `-addPL`)

The program can also simulate and print I16 tag (requires `-addI16`) and QS tags (requires `-addQS`).

You can also use the program to simulate unobserved inva

sites (requires `-explode`). This is useful for simulating genotype likelihoods for all sites in the input file, 



## Installation

```
git clone git@github.com:isinaltinkaya/vcfgl.git; cd vcfgl; make
```

## Usage


```
vcfgl

Usage: ./vcfgl -i <input> [options]


Options:
	-i/--input		Input file (required)
	-o/--output		Output file prefix (default:output)
	-O/--output-mode	Output mode (default:b)
				v	VCF file
				b	BCF file
				z	Compressed VCF file (vcf.gz)
				u	Uncompressed BCF file

	-e/--error-rate		Error rate (default:0.01)
	-d/--depth		Mean per-site read depth (default:1.0)
					Use `--depth inf` to set the simulated values to known true variables
	-df/--depths-file	File containing mean per-site read depth for each sample (conflicts with -d)
	--seed			Random seed used to initialize the random number generator
	--pos0			Are the input coordinates are 0-based? (default:0)
				If input cordinates are 0 based, use --pos0 1 to shift positions by +1

	-explode		Explode to unobserved sites in the input file (default:0=disabled)
	-addGP			Add GP field (default:0=disabled)
	-addPL			Add PL field (default:0=disabled)
	-addI16			Add I16 field (default:0=disabled)
	-addQS			Add QS field (default:0=disabled)

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
