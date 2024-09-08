# vcfgl - Miscellaneous Programs 

## fetchGl

fetchGl is one of the miscellaneous programs within the VCFGL project.

fetchGl is a command-line tool designed to fetch genotype likelihoods (GLs) corresponding to a specified genotype from a VCF file. It reads a VCF file, extracts the GLs for the specified genotype, and outputs them in a CSV format.

## Usage

```sh
fetchGl -i <input.vcf> -gt <genotype>
```

For example, 

```sh
./fetchGl -i input.vcf -gt AA
```

This command reads `input.vcf`, fetches the GLs for the genotype `AA`, and prints the results to the standard output.

## Installation

To compile the program, simply run `make` within the `misc/` directory. 

```sh
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl/misc;
make;
```

## Options

- `-i <input.vcf>`: Specifies the input VCF file.
- `-gt <genotype>`: Specifies the genotype to fetch. The genotype should be in the form of two alleles (e.g., `AA`, `AT`, `TT`).

## Output

The program outputs the GLs in CSV format with the following structure:

```
Position,GLs for samples
```

Each row represents the GLs for the specified genotype at a particular position in the VCF file.

## Example 


```sh
./misc/fetchGl -i test/testwd/test10.vcf -gt AC 2> stderr > stdout
```

```sh
$ cat stderr
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 1.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 2.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 3.
Requested genotype's first allele corresponds to allele 1 and second allele corresponds to allele 0 at position 4.
Requested genotype's first allele corresponds to allele 1 and second allele corresponds to allele 0 at position 5.
Requested genotype's first allele corresponds to allele 1 and second allele corresponds to allele 0 at position 6.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 7.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 8.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 9.
Requested genotype's first allele corresponds to allele 0 and second allele corresponds to allele 1 at position 10.
```

```sh
$ cat stdout
1,-0.596016,-0.297476
2,-1.190620,MISSING
3,-0.595338,-0.594828
4,-0.297687,0.000000
5,-0.893762,-0.297715
6,-0.298021,-0.595525
7,-0.297624,MISSING
8,-1.190520,MISSING
```

The MISSING value in the output indicates that the genotype likelihood (GL) value corresponding to the specified genotype for that particular sample and position is not available in the original VCF file. The `stderr` contains information about the allele indices corresponding to the requested genotype at each position. The `stdout` shows the GL values for each position and sample, with positions listed first, followed by the GL values.


### Comparison with Original VCF File

Here are the first two records from the original VCF file:

```
$ cat test/testwd/test10.vcf 
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	1	.	A	C,G,T,<*>	.	PASS	DP=3;QS=2,0,0,0,0;I16=2,1,0,0,49,801,0,0,60,1200,0,0,70,1650,0,0;AD=3,0,0,0,0;ADF=2,0,0,0,0;ADR=1,0,0,0,0	DP:GL:PL:GP:AD:ADF:ADR	2:0,-0.596016,-4.31724,-0.596016,-4.31724,-4.31724,-0.596016,-4.31724,-4.31724,-4.31724,-0.596016,-4.31724,-4.31724,-4.31724,-4.31724:0,6,43,6,43,43,6,43,43,43,6,43,43,43,43:1,0.253504,4.81682e-05,0.253504,4.81682e-05,4.81682e-05,0.253504,4.81682e-05,4.81682e-05,4.81682e-05,0.253504,4.81682e-05,4.81682e-05,4.81682e-05,4.81682e-05:2,0,0,0,0:2,0,0,0,0:0,0,0,0,0	1:0,-0.297476,-2.08531,-0.297476,-2.08531,-2.08531,-0.297476,-2.08531,-2.08531,-2.08531,-0.297476,-2.08531,-2.08531,-2.08531,-2.08531:0,3,21,3,21,21,3,21,21,21,3,21,21,21,21:1,0.504108,0.00821662,0.504108,0.00821662,0.00821662,0.504108,0.00821662,0.00821662,0.00821662,0.504108,0.00821662,0.00821662,0.00821662,0.00821662:1,0,0,0,0:0,0,0,0,0:1,0,0,0,0
ref1	2	.	A	C,G,T,<*>	.	PASS	DP=4;QS=1,0,0,0,0;I16=2,2,0,0,65,1059,0,0,80,1600,0,0,77,1639,0,0;AD=4,0,0,0,0;ADF=2,0,0,0,0;ADR=2,0,0,0,0	DP:GL:PL:GP:AD:ADF:ADR	4:0,-1.19062,-8.4554,-1.19062,-8.4554,-8.4554,-1.19062,-8.4554,-8.4554,-8.4554,-1.19062,-8.4554,-8.4554,-8.4554,-8.4554:0,12,85,12,85,85,12,85,85,85,12,85,85,85,85:1,0.0644729,3.50432e-09,0.0644729,3.50432e-09,3.50432e-09,0.0644729,3.50432e-09,3.50432e-09,3.50432e-09,0.0644729,3.50432e-09,3.50432e-09,3.50432e-09,3.50432e-09:4,0,0,0,0:2,0,0,0,0:2,0,0,0,0	0:.,.,.,.,.,.,.,.,.,.,.,.,.,.,.:.,.,.,.,.,.,.,.,.,.,.,.,.,.,.:.,.,.,.,.,.,.,.,.,.,.,.,.,.,.:0,0,0,0,0:0,0,0,0,0:0,0,0,0,0
```

From the above records,

- For position 1, `ind1` has GL `0,-0.596016,...` and `ind2` has GL  `0,-0.297476,...`. As we asked for the genotype `AC`, the GL values corresponding to this genotype are `-0.596016` for individual `ind1` and `-0.297476` for individual `ind2`, and reported as:

```
1,-0.596016,-0.297476
```

- For position 2, `ind1` has GL `0,-1.19062,...` and `ind2` has GL `.,.,.,.,.,.,.,.,.,.` indicating missing values, i.e. the individual `ind2` has no data at site. Genotype `AC` at this site is reported as:

```
2,-1.190620,MISSING
```

This explains why the output for position 2 includes `MISSING` for `ind2`.

Order of the genotypes [are given
by](https://github.com/samtools/htslib/blob/1187fa832998dd5fea9ea2a78bf6863a31c508f9/htslib/vcf.h#L1028-L1029):


```c
/** Conversion between alleles indexes to Number=G genotype index (assuming diploid, all 0-based) */
#define bcf_alleles2gt(a,b) ((a)>(b)?((a)*((a)+1)/2+(b)):((b)*((b)+1)/2+(a)))
```

Which can be summarised in a table as:

| a/b | 0   | 1   | 2   | 3   | 4   |
|-----|-----|-----|-----|-----|-----|
| **0** | 0   | 1   | 3   | 6   | 10  |
| **1** | 1   | 2   | 4   | 7   | 11  |
| **2** | 3   | 4   | 5   | 8   | 12  |
| **3** | 6   | 7   | 8   | 9   | 13  |
| **4** | 10  | 11  | 12  | 13  | 14  |

If REF is `A`, and ALTs are `C`,`G`,`T`,`<*>`:

| a/b  | A         | C         | G         | T         | <*>       |
|------|-----------|-----------|-----------|-----------|-----------|
| **A**    | AA (0)    | AC (1)    | AG (3)    | AT (6)    | A<*> (10) |
| **C**    | AC (1)    | CC (2)    | CG (4)    | CT (7)    | C<*> (11) |
| **G**    | AG (3)    | CG (4)    | GG (5)    | GT (8)    | G<*> (12) |
| **T**    | AT (6)    | CT (7)    | GT (8)    | TT (9)    | T<*> (13) |
| **<*>**  | A<*> (10) | C<*> (11) | G<*> (12) | T<*> (13) | <*><*> (14) |
