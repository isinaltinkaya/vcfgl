# Using vcfgl with msprime


We will use `msprime` to simulate a tree sequence with binary mutations and then use `vcfgl` to simulate genotype likelihoods and sequencing errors on the VCF file we generated. Let's start by simulating a tree sequence with binary mutations using `msprime`:

```python
import tskit
import msprime
ts = msprime.sim_ancestry(samples=3, ploidy=2, sequence_length=10, random_seed=2)
ts = msprime.sim_mutations(ts, rate=0.1, random_seed=2,model="binary")
with open("msprime_output.vcf","w") as vcffile:
    ts.write_vcf(vcffile,individual_names=["ind1","ind2","ind3"])
```

Let's examine the VCF file we simulated using `msprime`:

```bash
##fileformat=VCFv4.2
##source=tskit 0.5.6
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ind1    ind2    ind3
1       2       0       0       1       .       PASS    .       GT      0|0     0|0     0|0
1       3       1       0       1       .       PASS    .       GT      0|0     0|1     0|0
1       4       2       0       1       .       PASS    .       GT      1|0     1|1     0|0
1       6       3       0       1       .       PASS    .       GT      1|0     1|1     0|0
1       9       4       0       1       .       PASS    .       GT      0|0     1|0     0|0
```

Notice that the alleles in `REF` and `ALT` columns are binary (0 and 1). We will use `--source 0` with vcfgl to specify that the alleles are in binary format, which will assume that 0 is `A` and 1 is `C`.

Now, let's use `vcfgl` to simulate genotype likelihoods and sequencing errors on the VCF file we generated using `msprime`:

<details closed><summary><pre><code>$ ../vcfgl  -i msprime_output.vcf -O v -o sim_source0 --depth 10 --error-rate 0  --seed 42
        --depth 10 \
        --error-rate 0 \
        --source 0 \
        --seed 42
</code></pre></summary>
<pre>


vcfgl [version: v1.0.0-1648c26] [build: Apr  8 2024 16:35:03] [htslib: 1.15.1-20-g46c56fcc]


-&gt; Log file: sim_source0.arg
-&gt; Simulation output file: sim_source0.vcf

Command: vcfgl --verbose 0 --threads 1 --seed 42 --input msprime_output.vcf --source 0 --output sim_source0 --output-mode v --depth 10.000000 --error-rate 0.000000 --error-qs 0  --gl-model 2   --precise-gl 0   --adjust-qs 0  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 0 -printTruth 0  -printBasePickError 0 -printQsError 0 -printGlError 0 -printQScores 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 0 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0


[Program start] Mon Apr 08 16:38:44 2024


-&gt; Simulation finished successfully.

Summary:
	Number of samples: 3
	Number of contigs: 1
	Total number of sites simulated: 5
	Number of sites included in simulation output file: 5
	Number of sites skipped: 0

	Elapsed time (CPU): 0.000336 seconds
	Elapsed time (Real): 0.000000 seconds

</pre></pre></details>




<pre>$ cat sim_source0.vcf                                                                     
##fileformat=VCFv4.2
##FILTER=&lt;ID=PASS,Description=&quot;All filters passed&quot;&gt;
##source=tskit 0.5.6
##contig=&lt;ID=1,length=10&gt;
##fileDate=Mon Apr 08 16:38:44 2024
##source=vcfgl [version: v1.0.0-1648c26] [build: Apr  8 2024 16:35:03] [htslib: 1.15.1-20-g46c56fcc]
##source=Command: vcfgl --verbose 0 --threads 1 --seed 42 --input msprime_output.vcf --source 0 --output sim_source0 --output-mode v --depth 10.000000 --error-rate 0.000000 --error-qs 0  --gl-model 2   --precise-gl 0   --adjust-qs 0  -explode 0 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 -printPileup 0 -printTruth 0  -printBasePickError 0 -printQsError 0 -printGlError 0 -printQScores 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 1 -addInfoDP 0 -addFormatAD 0 -addInfoAD 0 -addFormatADF 0 -addInfoADF 0 -addFormatADR 0 -addInfoADR 0
##ALT=&lt;ID=*,Description=&quot;Symbolic alternate allele representing any possible alternative allele at this location&quot;&gt;
##FORMAT=&lt;ID=DP,Number=1,Type=Integer,Description=&quot;Simulated per-sample read depth&quot;&gt;
##FORMAT=&lt;ID=GL,Number=G,Type=Float,Description=&quot;Genotype likelihood in log10 likelihood ratio format&quot;&gt;
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2	ind3
1	2	0	A	&lt;*&gt;	.	PASS	.	DP:GL	11:0,-3.31133,-74.5483	11:0,-3.31133,-74.5483	8:0,-2.40824,-54.217
1	3	1	A	C,&lt;*&gt;	.	PASS	.	DP:GL	8:0,-2.40824,-54.217,-2.40824,-54.217,-54.217	11:-57.6828,0,-10.2429,-58.2848,-12.9522,-71.237	9:0,-2.70927,-60.9941,-2.70927,-60.9941,-60.9941
1	4	2	A	C,&lt;*&gt;	.	PASS	.	DP:GL	6:-11.7481,0,-25.3023,-12.9522,-25.9044,-38.8566	9:-60.9941,-2.70927,0,-60.9941,-2.70927,-60.9941	10:0,-3.0103,-67.7712,-3.0103,-67.7712,-67.7712
1	6	3	A	C,&lt;*&gt;	.	PASS	.	DP:GL	12:-30.2732,0,-43.8275,-32.3805,-45.3326,-77.7131	11:-74.5483,-3.31133,0,-74.5483,-3.31133,-74.5483	13:0,-3.91339,-88.1026,-3.91339,-88.1026,-88.1026
1	9	4	A	C,&lt;*&gt;	.	PASS	.	DP:GL	7:0,-2.10721,-47.4398,-2.10721,-47.4398,-47.4398	11:-37.3514,0,-30.5743,-38.8565,-32.3805,-71.237	16:0,-4.81648,-108.434,-4.81648,-108.434,-108.434
</pre>


N.B. If we did not use the binary model with msprime, the alleles would be in the regular VCF format (REF and ALT are one of `{A,C,G,T}`), and we would use `--source 1` with vcfgl.

---

For more information about the python code for obtaining VCF files please refer to <https://tskit.dev/tskit/docs/stable/python-api.html>.

---

