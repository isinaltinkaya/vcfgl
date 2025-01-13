# Simulating errors in the quality scores (`--error-qs`)

The `--error-qs` option in vcfgl allows the user to simulate errors in the quality scores. The option has three possible values:

```
   -eq, --error-qs [0]|1|2 __________ 0: Do not simulate errors in quality scores. Assumes all quality score assignments are correct
                                      1: Simulate site-specific errors in the probability of wrong base calls (requires: -bv FLOAT)
                                      2: Simulate the errors in the reported quality scores and genotype likelihoods (requires: -bv FLOAT)
```

## No errors in quality scores (`--error-qs 0`)

By default, vcfgl does not simulate errors in quality scores. This means the simulated quality scores are the values obtained from applying the phred scaling transformation to the specified base-calling error probability (`--error-rate`). 

Following is an example of simulating errors in quality scores with `--error-qs 0`, with an error rate of 0.4:

```
./vcfgl -i test/data/data2.vcf --depth 2 --error-rate 0.4 --error-qs 0 -addFormatAD 1 -printPileup 1 -s 42 -O v -printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1 -o error_qs0 > details_qs0.tsv
```

The options `-printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1` are used to print the error probabilities used in base picking process (`-printBasePickError 1`), the error probabilities used in quality score assignment (`-printQsError 1`), the error probabilities used in genotype likelihood calculations (`-printGlError 1`), and the quality scores (`-printQScores 1`), respectively. The results are printed to stdout, therefore we use `> details_qs0.tsv` to redirect the output to a file. The columns are: `type, sample_id, contig, site, read_index, value`


Let's check the output files:

```
$ cat details_qs0.tsv
base_pick_error_prob	NA	NA	NA	NA	0.400000
gl_error_prob	NA	NA	NA	NA	0.501187
qs	NA	NA	NA	NA	3
qs_error_prob	NA	NA	NA	NA	0.400000
```

We can see that the base picking error probability is 0.4, the genotype likelihood error probability is 0.501187, and the quality score error probability is 0.4. The quality score is 3, which corresponds to the Phred-33 encoded quality score (36) for the error probability of 0.4.


The resulting pileup file:


```
$ zcat error_qs0.pileup.gz
ref1	2	A	2	AC	$$	1	T	$
ref1	4	A	4	AAAC	$$$$	0	*	*
ref1	5	A	2	CC	$$	2	AT	$$
ref1	6	A	1	C	$	2	AA	$$
```

Notice that the quality score is " $ " for all reads, which corresponds to the Phred-33 encoded quality score (36) for the error probability of 0.4.


```math
-10 \times \log_{10}(0.4) = 3.97940
```

```math
3.9794  + 33 = 36.9794  \approx 36
```

```math
ASCII(36) = \$
```

```
$ cat error_qs0.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.4.1
##contig=<ID=ref1,length=10>
##fileDate=Sun Mar 31 14:21:27 2024
##source=vcfgl [version: v0.5-1b08819] [build: Mar 31 2024 14:21:01] [htslib: 1.15.1-20-g46c56fcc]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input ../test/data/data2.vcf --source 0 --output error_qs0 --output-mode v --depth 2.000000 --error-rate 0.400000 --error-qs 0  --gl-model 2   --precise-gl 0   --adjust-qs 0  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 0  -printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##ALT=<ID=*,Description="Symbolic alternate allele representing any possible alternative allele at this location">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	2	.	A	C,T,<*>	.	PASS	.	DP:GL:AD	2:-0.123909,0,-0.123909,-0.299484,-0.299484,-0.598968,-0.299484,-0.299484,-0.598968,-0.598968:1,1,0,0	1:-0.475059,-0.475059,-0.475059,-0.175575,-0.175575,0,-0.475059,-0.475059,-0.175575,-0.475059:0,0,1,0
ref1	4	.	A	C,<*>	.	PASS	.	DP:GL:AD	4:0,-0.22724,-0.950118,-0.526724,-1.12569,-1.42518:3,1,0	0:.,.,.,.,.,.:0,0,0
ref1	5	.	C	A,T,<*>	.	PASS	.	DP:GL:AD	2:0,-0.35115,-0.950118,-0.35115,-0.950118,-0.950118,-0.35115,-0.950118,-0.950118,-0.950118:2,0,0,0	2:-0.598968,-0.299484,-0.123909,-0.299484,0,-0.123909,-0.598968,-0.299484,-0.299484,-0.598968:0,1,1,0
ref1	6	.	A	C,<*>	.	PASS	.	DP:GL:AD	1:-0.475059,-0.175575,0,-0.475059,-0.175575,-0.475059:0,1,0	2:0,-0.35115,-0.950118,-0.35115,-0.950118,-0.950118:2,0,0
```




## Simulate site-specific errors in the probability of wrong base calls (`--error-qs 1`)

vcfgl can simulate the errors in the quality scores by simulating the errors in the actual base picking error probabilities for each site. This is done by sampling the probability of wrong base calls from a beta distribution, and using this probability for all reads from all individuals at site, making it site-specific. The error probabilities used in the quality score assignment and genotype likelihood calculations are not affected by this option. The error probabilities are sampled from a beta distribution with a mean equal to the one specified by the `--beta-variance` option. 

``` 
./vcfgl -i test/data/data2.vcf  --depth 2 --error-rate 0.4 --error-qs 1 --beta-variance 1e-1 -addFormatAD 1 -printPileup 1 -s 42 -O v -o error_qs1 -printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1 > details_qs1.tsv
```

Let's check the output files:

```
$ cat details_qs1.tsv
base_pick_error_prob	ind1	ref1	2	NA	0.621323
base_pick_error_prob	ind1	ref1	4	NA	0.960007
base_pick_error_prob	ind1	ref1	5	NA	0.525474
base_pick_error_prob	ind1	ref1	6	NA	0.002457
base_pick_error_prob	ind2	ref1	2	NA	0.621323
base_pick_error_prob	ind2	ref1	4	NA	0.960007
base_pick_error_prob	ind2	ref1	5	NA	0.525474
base_pick_error_prob	ind2	ref1	6	NA	0.002457
gl_error_prob	NA	NA	NA	NA	0.501187
qs	NA	NA	NA	NA	3
qs_error_prob	NA	NA	NA	NA	0.400000
```

We can see that the base picking error probabilities are different from the `--error-qs 0` example. The base picking error probability for the first sample at position 4 is 0.960007, meaning that the probability of the wrong base call is 0.960007. This is different from the `--error-qs 0` example, where the base picking error probability was 0.4 for all samples and sites.

```
$ zcat error_qs1.pileup.gz
ref1	2	A	2	AC	$$	1	T	$
ref1	4	A	4	CTGA	$$$$	0	*	*
ref1	5	A	2	CC	$$	2	AG	$$
ref1	6	A	1	C	$	2	AA	$$
```

Notice that the reported quality scores (`$`) are the same as the `--error-qs 0`, but some of the base calls are different. For example, compare the first read at position 4 of the first sample. The base call is now `C` instead of `A` (i.e. wrong base call).

No errors (`--error-qs 0`):

```
ref1	4	A	4	AAAC	$$$$	0	*	*
```

With errors (`--error-qs 1`):

```
ref1	4	A	4	CTGA	$$$$	0	*	*
```



Now, let's check the VCF:


```
$ cat error_qs1.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.4.1
##contig=<ID=ref1,length=10>
##fileDate=Sun Mar 31 14:21:27 2024
##source=vcfgl [version: v0.5-1b08819] [build: Mar 31 2024 14:21:01] [htslib: 1.15.1-20-g46c56fcc]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input ../test/data/data2.vcf --source 0 --output error_qs1 --output-mode v --depth 2.000000 --error-rate 0.400000 --error-qs 1 --beta-variance 1.000000e-01 --gl-model 2   --precise-gl 0   --adjust-qs 0  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 0  -printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##ALT=<ID=*,Description="Symbolic alternate allele representing any possible alternative allele at this location">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	2	.	A	C,T,<*>	.	PASS	.	DP:GL:AD	2:-0.123909,0,-0.123909,-0.299484,-0.299484,-0.598968,-0.299484,-0.299484,-0.598968,-0.598968:1,1,0,0	1:-0.475059,-0.475059,-0.475059,-0.175575,-0.175575,0,-0.475059,-0.475059,-0.175575,-0.475059:0,0,1,0
ref1	4	.	A	C,G,T,<*>	.	PASS	.	DP:GL:AD	4:-0.123909,0,-0.123909,-1.19209e-07,-5.96046e-08,-0.123909,0,-1.19209e-07,-5.96046e-08,-0.123909,-0.299484,-0.299484,-0.299484,-0.299484,-0.598968:1,1,1,1,0	0:.,.,.,.,.,.,.,.,.,.,.,.,.,.,.:0,0,0,0,0
ref1	5	.	C	A,G,<*>	.	PASS	.	DP:GL:AD	2:0,-0.35115,-0.950118,-0.35115,-0.950118,-0.950118,-0.35115,-0.950118,-0.950118,-0.950118:2,0,0,0	2:-0.598968,-0.299484,-0.123909,-0.299484,0,-0.123909,-0.598968,-0.299484,-0.299484,-0.598968:0,1,1,0
ref1	6	.	A	C,<*>	.	PASS	.	DP:GL:AD	1:-0.475059,-0.175575,0,-0.475059,-0.175575,-0.475059:0,1,0	2:0,-0.35115,-0.950118,-0.35115,-0.950118,-0.950118:2,0,0
```

With no errors example (`--error-qs 0`), we had the following for the contig ref1 at position 4:

```
ref1	4	.	A	C,<*>	.	PASS	.	DP:GL:AD	4:0,-0.22724,-0.950118,-0.526724,-1.12569,-1.42518:3,1,0	0:.,.,.,.,.,.:0,0,0
```


With errors (`--error-qs 1`) we have:

```
ref1	4	.	A	C,G,T,<*>	.	PASS	.	DP:GL:AD	4:-0.123909,0,-0.123909,-1.19209e-07,-5.96046e-08,-0.123909,0,-1.19209e-07,-5.96046e-08,-0.123909,-0.299484,-0.299484,-0.299484,-0.299484,-0.598968:1,1,1,1,0	0:.,.,.,.,.,.,.,.,.,.,.,.,.,.,.:0,0,0,0,0
```

The allelic read depth (AD) is different between the two VCF files. The AD for the first sample at position 4 is `3,1,0` for the no error example, meaning 3 reads with the reference allele (`A`), 1 read with the first alternate allele (`C`), and 0 reads with the second alternate allele (`<*>`). However, the AD is `1,1,1,1,0` for the `--error-qs 1` example, corresponding to 1 read with the reference allele (`A`), 1 read with the first alternate allele (`C`), 1 read with the second alternate allele (`G`), 1 read with the third alternate allele (`T`), and 0 reads with the fourth alternate allele (`<*>`). This is because the actual base calls were different between the two examples. 


## Simulate errors in quality scores through quality scores (`--error-qs 2`)

vcfgl can also simulate the errors in the quality scores by simulating the errors in the reported quality scores. This is done by sampling the error probability used in the reported quality scores from a beta distribution with a mean equal to the one specified by the `--beta-variance` option. For a more detailed description of this, please see vcfgl publication Supplementary Material, Algorithm 1.

``` 
vcfgl -i test/data/data2.vcf  --depth 2 --error-rate 0.4 --error-qs 2 --beta-variance 1e-1 -O v -addFormatAD 1 -printPileup 1 -o error_qs2 -s 42 -addFormatAD 1  -explode 1
```

Let's check the output files:

```
$ cat details_qs2.tsv
base_pick_error_prob	NA	NA	NA	NA	0.400000
gl_error_prob	ind1	ref1	2	0	0.630957
gl_error_prob	ind1	ref1	2	1	1.000000
gl_error_prob	ind1	ref1	4	0	0.002512
gl_error_prob	ind1	ref1	4	1	1.000000
gl_error_prob	ind1	ref1	4	2	0.794328
gl_error_prob	ind1	ref1	4	3	1.000000
gl_error_prob	ind1	ref1	5	0	0.050119
gl_error_prob	ind1	ref1	5	1	0.501187
gl_error_prob	ind1	ref1	6	0	0.063096
gl_error_prob	ind2	ref1	2	0	0.630957
gl_error_prob	ind2	ref1	5	0	0.794328
gl_error_prob	ind2	ref1	5	1	0.039811
gl_error_prob	ind2	ref1	6	0	0.199526
gl_error_prob	ind2	ref1	6	1	0.398107
qs	ind1	ref1	2	0	2
qs	ind1	ref1	2	1	0
qs	ind1	ref1	4	0	26
qs	ind1	ref1	4	1	0
qs	ind1	ref1	4	2	1
qs	ind1	ref1	4	3	0
qs	ind1	ref1	5	0	13
qs	ind1	ref1	5	1	3
qs	ind1	ref1	6	0	12
qs	ind2	ref1	2	0	2
qs	ind2	ref1	5	0	1
qs	ind2	ref1	5	1	14
qs	ind2	ref1	6	0	7
qs	ind2	ref1	6	1	4
qs_error_prob	ind1	ref1	2	0	0.621323
qs_error_prob	ind1	ref1	2	1	0.960007
qs_error_prob	ind1	ref1	4	0	0.002457
qs_error_prob	ind1	ref1	4	1	0.988480
qs_error_prob	ind1	ref1	4	2	0.753975
qs_error_prob	ind1	ref1	4	3	0.975551
qs_error_prob	ind1	ref1	5	0	0.046570
qs_error_prob	ind1	ref1	5	1	0.401989
qs_error_prob	ind1	ref1	6	0	0.059966
qs_error_prob	ind2	ref1	2	0	0.525474
qs_error_prob	ind2	ref1	5	0	0.700733
qs_error_prob	ind2	ref1	5	1	0.039285
qs_error_prob	ind2	ref1	6	0	0.160517
qs_error_prob	ind2	ref1	6	1	0.326249
```

We can see that only the base picking error probability is the same as the `--error-qs 0` example. The error probabilities used for assigning the quality scores and calculating the genotype likelihoods, as well as the actual quality scores, are different from both examples above.

```
$ zcat error_qs2.pileup.gz
ref1	2	A	2	AC	#!	1	T	#
ref1	4	A	4	AAAC	;!"!	0	*	*
ref1	5	A	2	CC	.$	2	AT	"/
ref1	6	A	1	C	-	2	AA	(%
```

Notice that the reported quality scores (`$`) are different from the `--error-qs 0`, but the actual base calls are the same. For example, compare the first read at position 4 of the first sample.

No errors (`--error-qs 0`):

```
ref1	4	A	4	AAAC	$$$$	0	*	*
```

With errors (`--error-qs 2`)

```
ref1	4	A	4	AAAC	;!"!	0	*	*
```
The reported quality score for the first read is now `;` instead of `$`, but the actual base call is still `A` (i.e. correct base call).


The VCF file:


```
$ cat error_qs2.vcf
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=tskit 0.4.1
##contig=<ID=ref1,length=10>
##fileDate=Sun Mar 31 14:21:27 2024
##source=vcfgl [version: v0.5-1b08819] [build: Mar 31 2024 14:21:01] [htslib: 1.15.1-20-g46c56fcc]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input ../test/data/data2.vcf --source 0 --output error_qs2 --output-mode v --depth 2.000000 --error-rate 0.400000 --error-qs 2 --beta-variance 1.000000e-01 --gl-model 2   --precise-gl 0   --adjust-qs 0  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 1 -printTruth 0  -printBasePickError 1 -printQsError 1 -printGlError 1 -printQScores 1 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 1 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##ALT=<ID=*,Description="Symbolic alternate allele representing any possible alternative allele at this location">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the REF and ALT alleles">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	2	.	A	C,T,<*>	.	PASS	.	DP:GL:AD	2:0,-0.406187,-inf,-0.105157,-0.545228,-0.244198,-0.105157,-0.545228,-0.244198,-0.244198:1,1,0,0	1:-0.244198,-0.244198,-0.244198,-0.105157,-0.105157,0,-0.244198,-0.244198,-0.105157,-0.244198:0,0,1,0
ref1	4	.	A	C,<*>	.	PASS	.	DP:GL:AD	4:-inf,-0.30103,-inf,0,-2.72397,-2.42294:3,1,0	0:.,.,.,.,.,.:0,0,0
ref1	5	.	C	A,T,<*>	.	PASS	.	DP:GL:AD	2:0,-0.469033,-2.22985,-0.469033,-2.22985,-2.22985,-0.469033,-2.22985,-2.22985,-2.22985:2,0,0,0	2:-1.85948,-1.91088,-1.96918,-0.295069,-0.346466,0,-1.85948,-1.91088,-0.295069,-1.85948:0,1,1,0
ref1	6	.	A	C,<*>	.	PASS	.	DP:GL:AD	1:-1.64882,-0.291389,0,-1.64882,-0.291389,-1.64882:0,1,0	2:0,-0.480868,-1.73711,-0.480868,-1.73711,-1.73711:2,0,0
```
