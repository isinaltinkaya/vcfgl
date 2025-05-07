# vcfgl - Miscellaneous Programs 

## setAlleles

setAlleles is a command-line tool designed to set the alleles in a VCF file. It reads a VCF file, modifies the alleles according to the specified reference and alternative alleles, and outputs the modified VCF file.

### Usage

```sh
Option descriptions:
     -s, --long-option TYPE [X] _____ Description
     -s                               Short option (if any)
         --long-option                Long option
                       TYPE           Type of the argument value, can be:
                                        - INT (integer)
                                        - STRING (string)
                                        - FILE (filename)
                                        - x|y|z (one of the listed values x, y or z)
                            [X]       Default argument value (if any)


Usage: setAlleles -i <input> [options]

    -h, --help _____________________  Print this help message and exit
    -v, --version __________________  Print version and build information and exit

Options:
    -i, --input FILE ________________ Input BCF file
    -a, --alleles FILE ______________ Alleles TSV file
    -O, --output-mode [b]|u|z|v _____ Output mode, can be:
                                        - b: Compressed BCF (.bcf)
                                        - u: uncompressed BCF (.bcf)
                                        - z: compressed VCF (.vcf.gz)
                                        - v: uncompressed VCF (.vcf)
    -o, --output STRING ['output'] __ Output filename prefix
```

### Input alleles TSV file

The input alleles TSV file should be formatted as a tab-separated values (TSV) file with two columns. The first column represents the reference allele, and the second column represents the alternative alleles. 

e.g. to set the reference allele to `A` and the alternative alleles to `T` and `G`, the TSV file should look like this:

```sh
A	T,G
```

- The alternative alleles can be one of the following: `A`, `C`, `G`, `T`, and `<*>` or `<NON_REF>`. 
- The program expects only one reference allele, at least one alternative allele and at most 4 alternative alleles. 
- Each line in the file corresponds to a site in the VCF file. The number of lines in the TSV file should be equal to the number of sites in the VCF file.


### Example

```sh
./setAlleles -i data/data2.vcf -a data/data2_alleles2.tsv -o output -O v
```

This command reads `data/data2.vcf`, modifies the alleles according to the specified reference and alternative alleles in `data/data2_alleles2.tsv`, and outputs the modified VCF file with the prefix `output` in uncompressed VCF format (i.e. `output.vcf`).

Example alleles TSV file:

```sh
$ cat data/data2_alleles2.tsv 
G	<*>,T
A	<*>
A	C
C	A
A	T,<*>,C,G
A	C
<*>	C,G,T,A
C	A
<*>	C,G,T,A
C	A,G,T
```

Example line from the input file:

```sh
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	1	.	A	C,T,G,<*>	.	PASS	DP=3;QS=0.5,0.5,1,0,0	DP:GL:PL:GP	2:0,-0.105667,0,-0.301034,-0.301034,-0.195367,-0.301034,-0.301034,-0.195367,-0.195367,-0.301034,-0.301034,-0.195367,-0.195367,-0.195367:0,1,0,3,3,2,3,3,2,2,3,3,2,2,2:0.104054,0.0815818,0.104054,0.0520268,0.0520268,0.0663581,0.0520268,0.0520268,0.0663581,0.0663581,0.0520268,0.0520268,0.0663581,0.0663581,0.0663581	1:-0.4,-0.4,-0.4,-0.301034,-0.301034,0,-0.4,-0.4,-0.301034,-0.4,-0.4,-0.4,-0.301034,-0.4,-0.4:4,4,4,3,3,0,4,4,3,4,4,4,3,4,4:0.0570268,0.0570268,0.0570268,0.0716218,0.0716218,0.143245,0.0570268,0.0570268,0.0716218,0.0570268,0.0570268,0.0570268,0.0716218,0.0570268,0.0570268
```

Corresponding line in the output file:

```sh
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
ref1	1	.	G	<*>,T	.	PASS	DP=3;QS=0,0,1	DP:GL:PL:GP	2:0,0,0,0,0,0:0,0,0,0,0,0:0.166667,0.166667,0.166667,0.166667,0.166667,0.166667	1:-0.4,-0.4,-0.4,-0.301034,-0.301034,0:4,4,4,3,3,0:0.12463,0.12463,0.12463,0.156527,0.156527,0.313057
```

Notice that, as described in the first line of the TSV file, the reference allele `A` has been replaced with `G`, and the alternative alleles `C`, `T`, `G`, `<*>` have been replaced with `<*>`, `T`, respectively. 


### Important notes

[!] The program only supports and updates the following INFO and FORMAT fields: `INFO/QS`, `FORMAT/GL`, `FORMAT/PL`, `FORMAT/GP`. 
[!] The program does not update any other fields in the VCF file.
[!] The program does not recalculate INFO/QS. It only updates the allele ordering in the INFO/QS field.
[!] The program does not recalculate DP values.

___

## fetchGl

fetchGl is a command-line tool designed to fetch genotype likelihoods (GLs) corresponding to a specified genotype from a VCF file. It reads a VCF file, extracts the GLs for the specified genotype, and outputs them in a CSV format.

### Usage

```sh
fetchGl -i <input.vcf> -gt <genotype>
```

For example, 

```sh
./fetchGl -i input.vcf -gt AA
```

This command reads `input.vcf`, fetches the GLs for the genotype `AA`, and prints the results to the standard output.

### Installation

To compile the program, simply run `make` within the `misc/` directory. 

```sh
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl/misc;
make;
```

### Options

- `-i <input.vcf>`: Specifies the input VCF file.
- `-gt <genotype>`: Specifies the genotype to fetch. The genotype should be in the form of two alleles (e.g., `AA`, `AT`, `TT`).

### Output

The program outputs the GLs in CSV format with the following structure:

```
Position,GLs for samples
```

Each row represents the GLs for the specified genotype at a particular position in the VCF file.

### Example 


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
