# vcfgl with SLiM

Let's perform a simple simulation with SLiM and then use `vcfgl` to simulate the sequencing data.

Example SLiM script `example.slim` contains the following code:

```
initialize() {
	setSeed(42);
	defineConstant("L", 20);
	initializeSLiMOptions(nucleotideBased=T);
	initializeAncestralNucleotides(randomNucleotides(L));
	initializeMutationTypeNuc("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(2.5e-5));
	initializeGenomicElement(g1, 0, L-1);
	initializeRecombinationRate(1e-8);
}

1 late() {
	defineConstant("simID",getSeed());
	sim.addSubpop("p1", 200);
	sim.addSubpop("p2", 200);
}

5000 late() {
	p1samples = p1.sampleIndividuals(5);
	p2samples = p2.sampleIndividuals(5);
	c(p1samples, p2samples).genomes.outputVCF(filePath="slim_out.vcf",simplifyNucleotides=T);
}
```

Run the SLiM script with the following command:

```
slim example.slim
```

The output file `slim_out.vcf` contains the following data:


```
##fileformat=VCFv4.2
##fileDate=20240403
##source=SLiM
##INFO=<ID=MID,Number=.,Type=Integer,Description="Mutation ID in SLiM">
##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">
##INFO=<ID=DOM,Number=.,Type=Float,Description="Dominance">
##INFO=<ID=PO,Number=.,Type=Integer,Description="Population of Origin">
##INFO=<ID=TO,Number=.,Type=Integer,Description="Tick of Origin">
##INFO=<ID=MT,Number=.,Type=Integer,Description="Mutation Type">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele Count">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=NONNUC,Number=0,Type=Flag,Description="Non-nucleotide-based">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1,URL=https://github.com/MesserLab/SLiM>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	i0	i1	i2	i3	i4	i5	i6	i7	i8	i9
1	3	.	T	C,G	1000	PASS	AC=2,8;DP=1000;AA=T	GT	0|0	0|0	0|0	0|0	0|0	1|2	1|2	2|2	2|2	2|2
1	7	.	T	G	1000	PASS	AC=10;DP=1000;AA=T	GT	0|0	0|0	0|0	0|0	0|0	1|1	1|1	1|1	1|1	1|1
1	8	.	A	C	1000	PASS	AC=9;DP=1000;AA=A	GT	0|0	0|0	0|0	0|0	0|0	1|1	0|1	1|1	1|1	1|1
1	10	.	T	A	1000	PASS	AC=2;DP=1000;AA=T	GT	0|1	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0
1	11	.	C	G,T	1000	PASS	AC=9,1;DP=1000;AA=C	GT	1|1	1|1	1|2	1|1	1|1	0|0	0|0	0|0	0|0	0|0
1	17	.	C	A,T	1000	PASS	AC=1,10;DP=1000;AA=C	GT	2|2	2|2	2|2	2|2	2|2	0|0	1|0	0|0	0|0	0|0
1	19	.	A	G	1000	PASS	AC=10;DP=1000;AA=A	GT	1|1	1|1	1|1	1|1	1|1	0|0	0|0	0|0	0|0	0|0
1	20	.	C	T	1000	PASS	AC=10;DP=1000;AA=C	GT	1|1	1|1	1|1	1|1	1|1	0|0	0|0	0|0	0|0	0|0
```

Now, let's use `vcfgl` to simulate the sequencing data. The following command will simulate the sequencing data with a depth of 1 and an error rate of 0.01:

```
vcfgl -i slim_out.vcf -d 1 -e 0.01 --source 1
```

