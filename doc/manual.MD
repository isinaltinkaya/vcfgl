# Errors in quality scores

You can use `--error-qs 1` to simulate errors in quality scores by simulating the site-specific errors in the probability of wrong base calls not accounted in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution-sampled error rates and the reported quality scores use the beta distribution mean (i.e. the error rate).

You can use `--error-qs 2` to simulate errors in quality scores by simulating the errors in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution mean (i.e. the error rate) and the reported quality scores use the beta distribution-sampled error rates.


<!-- 
## Example

```shell
./vcfgl -i test/t1.vcf -o test/t1_out -e 0.01 --seed 42 --depth 1
```

Input file: `test/t1.vcf`

```
##fileformat=VCFv4.2
##source=tskit 0.4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind
chr22	2	.	0	1	.	PASS	.	GT	0|0
chr22	10	.	0	1	.	PASS	.	GT	1|0
chr22	42	.	0	1	.	PASS	.	GT	0|1
chr22	98	.	0	1	.	PASS	.	GT	1|1
```

Output file: `test/t1_out.bcf`

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.4.1
##contig=<ID=chr22,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##fileDate=Tue Nov 21 09:56:50 2023
##source=vcfgl [version: v0.3.3 bb8df2d-dirty] [build: Nov 21 2023 09:56:39] [htslib: 1.15.1-20-g46c56fc]
##source=vcfgl --verbose 0 --threads 1 --input test/t1.vcf --output test/t1_out --output-mode b --depth 1.000000 --error-rate 0.010000 --error-qs 0 --beta-variance -1.000000e+00 --precise-gl 1 --seed 42 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -explode 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 0 -addFormatAD 0 -addFormatADF 0 -addFormatADR 0 -addInfoDP 0 -addInfoAD 0 -addInfoADF 0 -addInfoADR 0
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##bcftools_viewVersion=1.18+htslib-1.18
##bcftools_viewCommand=view test/t1_out.bcf; Date=Tue Nov 21 09:56:53 2023
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind
chr22	2	.	A	C,G,T	.	PASS	.	GT:GL	0|0:0,-0.29957,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276,-2.47276
chr22	10	.	A	C,G,T	.	PASS	.	GT:GL	1|0:.,.,.,.,.,.,.,.,.,.
chr22	42	.	A	C,G,T	.	PASS	.	GT:GL	0|1:-1.87362,0,-1.87362,-2.17319,-2.17319,-4.34637,-2.17319,-2.17319,-4.34637,-4.34637
chr22	98	.	C	A,G,T	.	PASS	.	GT:GL	0|0:0,-0.29957,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276,-2.47276
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
``` -->

--- 

# Ancestral and derived alleles

When `--source 0` is used, `vcfgl` treats the allelic state `0` in input file as `A`, and `1` as `C`. Therefore a genotype of `00` corresponds to `AA`, `01` to `AC`, `10` to `CA`, and `11` to `CC`.

Binary representation | Allele
-- | --
0 | A
1 | C

---

# Sampling strands

The following options require strand information:
* `-addI16 1`
* `-addFormatADF 1`
* `-addFormatADR 1`
* `-addInfoADF 1`
* `-addInfoADR 1`

These options will trigger the random sampling of strands. The strands are sampled as forward and reverse strands from a binomial distribution with equal probability.

---

# Unobserved allele notation

Unobserved alleles are defined as the alleles that are not observed at a site, and is represented by `<*>` in BCFtools and `NON_REF` in GATK. In vcfgl, there are 6 options for handling unobserved alleles:

* `-doUnobserved 0`: Trim the unobserved alleles
* `-doUnobserved 1`: Use symbolic alternate allele notation `<*>` for unobserved alleles. This notation is often used by BCFtools.
* `-doUnobserved 2`: Use `<NON_REF>` notation for unobserved alleles. This notation is often used by GATK.
* `-doUnobserved 3`: Enumerate the unobserved alleles as `A`, `C`, `G`, and `T`.
* `-doUnobserved 4`: Use symbolic alternate allele notation `<*>` for unobserved alleles and enumerate the unobserved alleles as `A`, `C`, `G`, and `T`.
* `-doUnobserved 5`: Use `<NON_REF>` notation for unobserved alleles and enumerate the unobserved alleles as `A`, `C`, `G`, and `T`.

For example, if we observe 10 reads for `A` and 4 reads for `C` at a site, the following table shows the output for different `-doUnobserved` options:

Option | REF field | ALT field
-- | -- | --
`-doUnobserved 0` | `A` | `C`
`-doUnobserved 1` | `A` | `C,<*>`
`-doUnobserved 2` | `A` | `C,<NON_REF>`
`-doUnobserved 3` | `A` | `C,G,T`
`-doUnobserved 4` | `A` | `C,G,T,<*>`
`-doUnobserved 5` | `A` | `C,G,T,<NON_REF>` -->