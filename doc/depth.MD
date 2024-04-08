
# Simulating read depths


---

- `-d, --depth FLOAT|'inf'`: Specify the mean per-site read depth. If 'inf' is provided, the read depths are assumed to be infinite and no read depth simulation is performed. Instead, the genotype likelihoods are calculated based on the true genotypes in the VCF file assuming infinite read depths.
- `-df, --depths-file FILE`: Specify per-sample mean per-site read depth values (one value per line)
---

vcfgl simulates read depths for each site and each sample in the VCF file. The read depths are simulated from a Poisson distribution with mean equal to the given depth value. The procedure is as follows:

1. For each site; for each haplotype; sample number of reads from Poisson distribution with mean equal to given depth value

2. Simulating error

(2.1) Probability of error: e

(2.2) Probability of sampling an incorrect allele due to error: error rate × 1/3

(2.3) Probability of sampling correct allele: 1 − error rate

4. Normalize the genotype likelihoods by substracting the maximum genotype likelihood value
observed at site.


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
```shell
./vcfgl -i input.vcf --depths-file depths.txt
```


## Set the values to known true variables using `--depth inf`

This is useful for comparing the genotype calling-based methods with genotype likelihood-based methods. It will set the values in GL, GP, and PL tags to the best possible values corresponding to the input genotypes. Therefore, instead of using `--depth 100`, you can simulate the best possible values quickly.
