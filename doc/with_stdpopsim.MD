# Using vcfgl with stdpopsim

We will use `stdpopsim` to simulate a tree sequence and then use `vcfgl` to simulate genotype likelihoods and sequencing errors on the VCF file we generated.

```python
import stdpopsim #v0.2.1

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("OutOfAfrica_3G09")
contig = species.get_contig("chr22",mutation_rate=model.mutation_rate)
samples = {"YRI": 2, "CHB": 2, "CEU": 0}
engine = stdpopsim.get_engine("msprime")
ts = engine.simulate(model, contig, samples)
with open("output.trees", "w") as treefile:
    ts.dump(treefile)
```

After obtaining the tree sequence, we can get VCF file using one of the following methods:

1. Using the tree sequence object `ts` with Python tskit library:

```python
with open("stdpopsim_output.vcf", "w") as vcf_file:
    ts.write_vcf(vcf_file)
```

2. Using tskit command line tool:

```bash
tskit vcf output.trees > stdpopsim_output.vcf
```

```bash
$ head stdpopsim_output.vcf
##fileformat=VCFv4.2
##source=tskit 0.5.6
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=50818468>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  tsk_0   tsk_1   tsk_2   tsk_3
1       53      0       C       G       .       PASS    .       GT      0|0     0|0     1|0     1|0
1       558     1       T       A       .       PASS    .       GT      1|0     0|1     0|0     0|0
1       1093    2       A       C       .       PASS    .       GT      0|1     0|0     0|0     0|0
1       1319    3       C       A       .       PASS    .       GT      1|0     0|1     0|0     0|0
```

Now, let's use `vcfgl` to simulate genotype likelihoods and sequencing errors on the VCF file we generated.

<details closed><summary><pre><code>$ vcfgl -i stdpopsim_output.vcf -O b -o sim \
        --depth 1 \
        --error-rate 0.002 \
        --source 1 \
        --seed 42 \
        -printTruth 1 \
        -addFormatAD 1  \
        -printPileup 1 
</code></pre></summary> 
<pre>
<pre>
vcfgl [version: v0.4.3-9dd8010] [build: Feb  3 2024 11:08:53] [htslib: 1.16]

Command: vcfgl --verbose 0 --threads 1 --seed 42 --input stdpopsim_output.vcf --source 1 --output output --output-mode b --depth 1.000000 --error-rate 0.002000 --error-qs 0 --gl-model 2 --platform 0 --precise-gl 0 --i16-mapq 20 -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 1 -printBasePickError 0 -printQsError 0 -printGlError 0 -printQScores 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0

Wed Apr 3 16:06:10 2024

-&gt; Log file: output.arg
-&gt; Simulation output file: output.bcf
-&gt; True genotypes output file: output.truth.bcf
-&gt; Pileup output file: output.pileup.gz

-&gt; Simulation finished successfully.

Summary:
Number of samples: 4
Number of contigs: 1
Total number of sites simulated: 103578
Number of sites included in simulation output file: 103578
Number of sites skipped: 0

</pre>
</pre> </details>

- Let's examine the command we used:

  - `-i` specifies the input VCF file. In this example, we use the VCF file `stdpopsim_output.vcf` that we previously generated as input to `vcfgl`.

  - `-O` specifies the output mode. In this example, `b` is used to output the results in BCF format (`.bcf` file).

  - `-o` specifies the output filename prefix. In this example, the output files will be named `output.bcf`

  - `--depth` specifies the average read depth per site. This is used as the mean of the Poisson distribution to simulate the read depth for each individual at each site.

  - `--error-rate` specifies the sequencing error rate.

  - `--source` specifies the format of the alleles in the input VCF file.

    - If the alleles are in binary format (REF=0, ALT=1), then `--source 0` should be used. This is typically the output format if you use the msprime BinaryMutationModel.
    - If the alleles are in regular VCF format (REF and ALT are one of `{A,C,G,T}`), then `--source 1` should be used. In this example, we use `--source 1` since we are using the regular VCF format.

- `--seed` specifies the seed for the random number generator. This is important for ensuring reproducibility of the simulation results. We use 42 in this example as it is the answer to the ultimate question of life, the universe, and everything.

- `-printTruth 1` specifies whether to output the true genotypes in the output file. In this example, we set it to 1 to output the true genotypes, which will be stored in the `output.truth.bcf` file.

```bash
$ bcftools view output.truth.bcf |head -15
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.5.6
##contig=<ID=1,length=50818468>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##fileDate=Wed Apr  3 16:06:10 2024
##source=vcfgl [version: v0.4.3-9dd8010] [build: Feb  3 2024 11:08:53] [htslib: 1.16]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input stdpopsim_output.vcf --source 1 --output output --output-mode b --depth 1.000000 --error-rate 0.002000 --error-qs 0  --gl-model 2 --platform 0 --precise-gl 0 --i16-mapq 20  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 1  -printBasePickError 0 -printQsError 0 -printGlError 0 -printQScores 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##bcftools_viewVersion=1.18-6-gc7cbe0b3+htslib-1.18-7-gac70212
##bcftools_viewCommand=view output.truth.bcf; Date=Wed Apr  3 16:13:28 2024
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tsk_0	tsk_1	tsk_2	tsk_3
1	53	    0	C	G	.	PASS	.	GT	0|0	0|0	1|0	1|0
1	558	    1	T	A	.	PASS	.	GT	1|0	0|1	0|0	0|0
1	1093	2	A	C	.	PASS	.	GT	0|1	0|0	0|0	0|0
1	1319	3	C	A	.	PASS	.	GT	1|0	0|1	0|0	0|0
```

Let's examine the simulated data in the `output.bcf` file:

```bash
$ bcftools view output.bcf |head -18
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.5.6
##contig=<ID=1,length=50818468>
##fileDate=Wed Apr  3 16:06:10 2024
##source=vcfgl [version: v0.4.3-9dd8010] [build: Feb  3 2024 11:08:53] [htslib: 1.16]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input stdpopsim_output.vcf --source 1 --output sim --output-mode b --depth 1.000000 --error-rate 0.002000 --error-qs 0  --gl-model 2 --platform 0 --precise-gl 0 --i16-mapq 20  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 1  -printBasePickError 0 -printQsError 0 -printGlError 0 -printQScores 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##ALT=<ID=*,Description="Symbolic alternate allele representing any possible alternative allele at this location">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the REF and ALT alleles">
##bcftools_viewVersion=1.18-6-gc7cbe0b3+htslib-1.18-7-gac70212
##bcftools_viewCommand=view output.bcf; Date=Wed Apr  3 16:16:08 2024
#CHROM POS  ID	REF ALT QUAL	FILTER	INFO	FORMAT	tsk_0	tsk_1	tsk_2	tsk_3
1	53	    0	C	G,<*>   	.	PASS	.	DP:GL:AD	1:0,-0.300741,-3.17625,-0.300741,-3.17625,-3.17625:1,0,0	0:.,.,.,.,.,.:0,0,0	1:0,-0.300741,-3.17625,-0.300741,-3.17625,-3.17625:1,0,0	2:-2.57477,0,-2.57477,-2.87551,-2.87551,-5.75103:1,1,0
1	558	    1	T	A,<*>	    .	PASS	.	DP:GL:AD	1:0,-0.300741,-3.17625,-0.300741,-3.17625,-3.17625:1,0,0	3:-5.45028,0,-2.27403,-5.75103,-2.87551,-8.62654:1,2,0	1:0,-0.300741,-3.17625,-0.300741,-3.17625,-3.17625:1,0,0	1:0,-0.300741,-3.17625,-0.300741,-3.17625,-3.17625:1,0,0
1	1093	2	A	<*>	        .	PASS	.	DP:GL:AD	0:.,.,.:0,0	0:.,.,.:0,0	0:.,.,.:0,0	1:0,-0.300741,-3.17625:1,0
1	1319	3	A	C,<*>	    .	PASS	.	DP:GL:AD	0:.,.,.,.,.,.:0,0,0	2:0,-0.601481,-6.35251,-0.601481,-6.35251,-6.35251:2,0,0	2:-6.35251,-0.601481,0,-6.35251,-0.601481,-6.35251:0,2,0	0:.,.,.,.,.,.:0,0,0
```

---

For more information about the python code for obtaining VCF files please refer to <https://tskit.dev/tskit/docs/stable/python-api.html>.

---
