# vcfgl
## Small program to simulate genotype likelihoods from VCF GT tags


```
./vcfgl -in {INPUT_VCF} -out {OUTPUT_FILE_PREFIX} -err {ERROR_RATE} -depth {DEPTH} -seed {SEED}
```

Arguments:

* **INPUT_VCF**

 
Input VCF file obtained from `tskit.write_vcf`.

After simulating a Tree Sequence, save your tree sequence as a VCF file:
```
with open("simulated.vcf", "w") as vcf_file:
    tree_sequence.write_vcf(vcf_file, ploidy=2)
```

For more information please refer to https://tskit.dev/tskit/docs/stable/python-api.html.


This program uses the GT tags in `simulated.vcf` to simulate genotype likelihoods.


* **OUTPUT_FILE_PREFIX**

File prefix for the output VCF and arguments file the program will generate. e.g. Using `-out myout` will generate the files `myout.vcf` and `myout.arg`.

* **ERROR_RATE**

Error rate to be used in simulation. For more details, see [Genotype likelihood calculation](#genotype-likelihood-calculation).

* **DEPTH**

Mean per site read depth to be used in simulations.

* **SEED**

Random seed used to initialize the random number generator.

## Genotype likelihood calculation

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

## Usage

```
git clone git@github.com:isinaltinkaya/vcf-gl.git; cd vcf-gl; make
```

