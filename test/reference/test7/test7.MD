# Test 7: gVCF test with 1 sample

```
$./vcfgl -i test/data/data1.vcf --verbose 0 --threads 1 --seed 42 --output-mode v
--depth 1 --error-rate 0.2 --error-qs 0  -GL 1 --gl1-theta 0.83 --platform 0
--precise-gl 0 -explode 1 --rm-invar-sites 0 --rm-empty-sites 1 -doUnobserved 2
-doPileup 1 -doTruth 1 -addGL 1 -addGP 1 -addPL 1 -addI16 0 -addQS 1
-addFormatDP 1 -addInfoDP 1 -addFormatAD 1 -addInfoAD 1 -addFormatADF 1
-addInfoADF 1 -addFormatADR 1 -addInfoADR 1 -dogvcf 1 --gvcf-dps 1 -o test7

$ cat test7.vcf | bcftools query -f "[%SAMPLE\t%CHROM,%POS,%REF,%ALT\tEND=%END\tMIN_DP=%MIN_DP\tPL=%PL\tDP=%DP\n]"
ind1	ref1,1,A,<NON_REF>	END=3	MIN_DP=1	PL=0,3,6	DP=1
ind1	ref1,4,A,C,<NON_REF>	END=4	MIN_DP=.	PL=1,0,1,4,4,4	DP=2
ind1	ref1,5,A,C,T,<NON_REF>	END=5	MIN_DP=.	PL=0,1,0,1,1,0,3,3,3,2	DP=3
ind1	ref1,7,A,<NON_REF>	END=8	MIN_DP=1	PL=0,3,6	DP=1
ind1	ref1,10,A,<NON_REF>	END=10	MIN_DP=1	PL=0,3,6	DP=1
```

Bcftools reference:

```
$ ../bcftools/bcftools mpileup test/data/data1.sam  --fasta-ref test/data/data1_ref.fa  -g 1 | bcftools query -f "[%SAMPLE\t%CHROM,%POS,%REF,%ALT\tEND=%END\tMIN_DP=%MIN_DP\tPL=%PL\tDP=%DP\n]"
[warning] The -a DP option is required with --gvcf, switching on.
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 250
test/data/data1.sam	ref1,1,A,<*>	END=3	MIN_DP=1	PL=0,3,6	DP=1
test/data/data1.sam	ref1,4,A,C,<*>	END=4	MIN_DP=.	PL=1,0,1,4,4,4	DP=2
test/data/data1.sam	ref1,5,A,T,C,<*>	END=5	MIN_DP=.	PL=0,1,0,1,1,0,3,3,3,2	DP=3
test/data/data1.sam	ref1,7,A,<*>	END=8	MIN_DP=1	PL=0,3,6	DP=1
test/data/data1.sam	ref1,10,A,<*>	END=10	MIN_DP=1	PL=0,3,6	DP=1
```

Can also be found at:

```
$ bcftools mpileup test/data/data1.sam  --fasta-ref test/data/data1_ref.fa  -g 1 | bcftools query -f "[%SAMPLE\t%CHROM,%POS,%REF,%ALT\tEND=%END\tMIN_DP=%MIN_DP\tPL=%PL\tDP=%DP\n]" > test/reference/test7/test7_gvcf_g1_query.tsv
```
