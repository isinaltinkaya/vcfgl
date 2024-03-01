<!--
TODO:
- [ ] Move examples etc (commented out stuff at the end of this file) to documentation
- [ ] Add issue template etc
- [ ] Add tutorials links in readme

-->

<a name="readme-top"></a>

<h3 align="center">vcfgl</h3>

<p align="center">
<a href="https://github.com/isinaltinkaya/vcfgl/actions/workflows/test.yml"><img src="https://github.com/isinaltinkaya/vcfgl/actions/workflows/test.yml/badge.svg" /></a> 
<a href="https://github.com/isinaltinkaya/vcfgl/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GNU%20GPLv3.0-purple.svg"/></a>
</p>
<!-- ![vcfgl](https://img.shields.io/badge/version-v0.3.3-brightgreen.svg)  -->

  <p align="center">
    Genotype likelihood simulator for VCF/BCF files
    <br />
    <!-- <a href="https://github.com/isinaltinkaya/vcfgl"><strong>Quickstart»</strong></a> -->
    <br />
    <br />
    <a href="https://github.com/isinaltinkaya/vcfgl">Installation</a>
    ·
    <a href="https://github.com/isinaltinkaya/vcfgl">Tutorial</a>
    ·
    <a href="https://github.com/isinaltinkaya/vcfgl/issues">Report Bug</a>
    ·
    <a href="https://github.com/isinaltinkaya/vcfgl/issues">Request Feature</a>
  </p>
</div>

<details open>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li>
      <a href="#installation">Installation</a>
      <ul>
        <!-- <li><a href="#prerequisites">Prerequisites</a></li> -->
        <!-- <li><a href="#method-1-using-htslib-submodule-make">Method 1: Using HTSlib submodule (`make`)</a></li> -->
        <!-- <li><a href="#method-2-using-systemwide-htslib-installation-make-htssrcsystemwide">Method 2: Using systemwide HTSlib installation (`make HTSSRC=systemwide`)</a></li> -->
        <!-- <li><a href="#method-3-using-specified-htslib-path-make-htssrcpathtohtslib">Method 3: Using specified HTSlib path (`make HTSSRC=path/to/htslib`)</a></li> -->
      </ul>
    </li>
    <li><a href="#usage">Usage</a>
    <ol>
      <li><a href="#general-options">General Options</a></li>
      <li><a href="#inputoutput">Input/Output</a></li>
      <li><a href="#simulation-parameters">Simulation Parameters</a></li>
      <li><a href="#output-vcfbcf-tags">Output VCF/BCF Tags</a></li>
    </ol>
    </li>
      <li><a href="#tutorials">Tutorials</a>
      <ol>
        <li><a href="#doGVCF">Simulating gVCF files</a></li>
      </ol>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#how-to-cite">How to cite</a></li>
  </ol>
</details>

## Overview

**vcfgl** is a lightweight command-line program for simulating VCF/BCF and gVCF files. It allows you to simulate sequencing data with various parameters, such as read depth, base-calling error rates, quality score errors, and genotype likelihood models.

## Installation

<!-- ### Prerequisites

HTSlib -->

You can install **vcfgl** using one of the following methods:

### &rarr; Method 1: Using HTSlib submodule

This method uses the htslib submodule in the repository.

```shell
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl;
make;
```

### &rarr; Method 2: Using systemwide HTSlib installation

This method assumes you have a systemwide htslib installation.

```shell
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl;
make HTSSRC="systemwide";
```

### &rarr; Method 3: Using specified HTSlib path

This method allows you to specify the path to your htslib installation.

```shell
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl;
make HTSSRC=path/to/htslib;
```

For detailed compilation instructions, in vcfgl directory, run:

<details><summary> <code> make help </code></summary>

<pre>$ make help

----------------------------------------
 Program: vcfgl
 Version: v0.4-bb8df2d
 License: GNU GPLv3.0
----------------------------------------

 Usage:
   make [target] [FLAG=value...]

 Targets:
   help    - Print this help message
   dev     - Compile in developer/debug mode (activates flags: -g -Wall -O0)
   clean   - Clean up the directory
   test    - Run unit tests

 Flags:
   HTSSRC  - Specifies the source of HTSlib.
       Values:
       (empty)          - Use the HTSlib submodule [default]
       systemwide       - Use the systemwide HTSlib installation
       /path/to/htslib  - Use the HTSlib installation at /path/to/htslib

 Examples:
   make                         - Compile in release mode using HTSlib submodule
   make HTSSRC=/path/to/htslib  - Compile in release mode using /path/to/htslib
   make dev HTSSRC=systemwide   - Compile in developer mode using the systemwide HTSlib installation

 Note: If no values are provided for HTSSRC, CXX, CXXFLAGS, or LIBS, defaults will be used.

</pre>
</details>

# Usage

You can access the command-line help page using `vcfgl --help` or `vcfgl -h`:

<details open>
<summary> <code> vcfgl --help </code></summary>

<pre>$ ./vcfgl -h

Usage: vcfgl -i &lt;input&gt; [options]

Option descriptions:
     -s, --long-option TYPE [X] _____ Description
     -s                               Short option (if any)
         --long-option                Long option
                       TYPE           Type of the argument value, can be:
                                        - INT (integer)
                                        - INT+ (additive integer: sum values to use together
                                        - FLOAT (floating point number)
                                        - STRING (string)
                                        - FILE (filename)
                                        - x|y|z (one of the listed values x, y or z)
                            [X]       Default argument value (if any)
                                _____ Connector to the option description for better readability


General options:
    -V, --verbose INT [0] ___________ Verbosity level
    -@, --threads INT [1] ___________ Number of threads
    -s, --seed INT [time] ___________ Random seed for initializing the random number generator

Input/Output:
    -i, --input FILE ________________ Input VCF/BCF file
        --source [0]|1 ______________ 0: Input REF/ALT alleles are in binary format (REF=0, ALT=1; typically outputted from msprime BinaryMutationModel)
                                      1: Input REF/ALT alleles are in VCF format (REF=i, ALT=j(,k..); i, j and k from {A,C,G,T}; i.e. the regular VCF format)
    -o, --output STRING [&apos;output&apos;] __ Output filename prefix
    -O, --output-mode [b]|u|z|v _____ b: Compressed BCF (.bcf), u: uncompressed BCF (.bcf), z: compressed VCF (.vcf.gz), v: uncompressed VCF (.vcf)

Simulation parameters:
    -d, --depth FLOAT|&apos;inf&apos; _________ Mean per-site read depth
                                      (&apos;inf&apos;) Simulate true values (requires: -addFormatDP 0 -addInfoDP 0)
   -df, --depths-file FILE __________ File containing mean per-site read depth values for each sample. One value per line.
    -e, --error-rate FLOAT __________ Base-calling error probability
   -eq, --error-qs [0]|1|2 __________ 0: Do not simulate errors in quality scores. Assumes all quality score assignments are correct
                                      1: Simulate site-specific errors in the probability of wrong base calls (requires: -bv FLOAT)
                                      2: Simulate the errors in the reported quality scores and genotype likelihoods (requires: -bv FLOAT)
   -bv, --beta-variance FLOAT _______ Designated variance for the beta distribution
   -GL, --gl-model 1|[2] ____________ Genotype likelihood model to be used in simulation
                                      1: Genotype likelihood model with correlated errors (a.k.a. Li model, samtools model, angsd -GL 1)
                                      2: Canonical genotype likelihood model with independent errors (a.k.a. McKenna model, GATK model, angsd -GL 2)
         --gl1-theta FLOAT [0.83] ___ Theta parameter for the genotype likelihood model 1 (requires: -GL 1)
         --platform [0]|1 ___________ 0: Do not use platform specification
                                      1: NovaSeq 6000 (RTA3), qualities are binned into 4 values: 2, 12, 23 and 37
         --precise-gl [0]|1 _________ 0: Use the discrete phred-scaled error probabilities in the genotype likelihood calculation
                                      1: Use precise error probabilities in the genotype likelihood calculation (requires: -GL 2)
         --i16-mapq INT [20] ________ Mapping quality score for I16 tag (requires: -addI16 1)
         --gvcf-dps INT(,INT..) _____ Minimum per-sample read depth range(s) for constructing gVCF blocks (requires: -doGVCF 1)
                                      Example: `--gvcf-dps 5,10,20` will group invariable sites into three types of gVCF blocks: [5,10), [10,20) and [20,inf)
                                      Sites with minimum depth &lt; 5 will be printed as regular VCF records.
         --adjust-qs INT+ [3] _______ 0: Do not adjust quality scores
                                      1: Use adjusted quality scores in genotype likelihoods (requires: --precise-gl 0)
                                      2: Use adjusted quality scores in calculating the quality score sum (QS) tag (requires: -addQS 1)
                                      4: Use adjusted quality scores in pileup output (requires: --printPileup 1)
                                      8: Use adjusted quality scores in --printQScores output (requires: --printQScores 1)
                                      16: Use adjusted quality scores in --printGlError output (requires: --printGlError 1)
         --adjust-by FLOAT [0.499] __ Adjustment value for quality scores (requires: --adjust-qs &gt; 0)
         -explode [0]|1 _____________ 1: Explode to sites that are not in input file.
                                      Useful for simulating invariable sites when the input file only contains variable sites.
                                      Sets all genotypes in exploded sites to homozygous reference.
         --rm-invar-sites INT+ [0]___ 0: Do not remove invariable sites
                                      1: Remove sites where all individuals&apos; true genotypes in the input file are homozygous reference
                                      2: Remove sites where all individuals&apos; true genotypes in the input file are homozygous alternative
                                      4: Remove sites where the all simulated reads among all individuals are the same base
                                      Example: &apos;--rm-invar-sites 3&apos; (1+2) will do both 1 and 2 (i.e. remove all homozygous sites)
         --rm-empty-sites [0]|1 _____ 0: Do not remove empty sites
                                      1: Remove empty sites (i.e. sites where no reads were simulated)
         -doUnobserved INT [1] ______ 0: Trim unobserved alleles. Only alleles that are observed will be listed in REF and ALT fields
                                      1: Use &apos;&lt;*&gt;&apos; notation to represent unobserved alleles
                                      2: Use &apos;&lt;NON_REF&gt;&apos; notation to represent unobserved alleles (a.k.a. GATK notation)
                                      3: Explode unobserved bases from {A,C,G,T} list
                                      4: Use &apos;&lt;*&gt;&apos; notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list
                                      5: Use &apos;&lt;NON_REF&gt;&apos; notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list
         -doGVCF [0]|1 ______________ 0: Disabled, 1: Output in gVCF format (requires: --rm-invar-sites 0, -doUnobserved 2, -addPL 1 and --gvcf-dps INT)
         -printPileup [0]|1 _________ 0: Disabled, 1: Also output in pileup format (&lt;output_prefix&gt;.pileup.gz)
         -printTruth [0]|1 __________ 0: Disabled, 1: Also output the VCF file containing the true genotypes (named &lt;output_prefix&gt;.truth.vcf)
         -printBasePickError [0]|1 __ 0: Disabled, 1: Print the base picking error probability to stdout.
                                      If --error-qs 1 is used, writes per-read base picking error probabilities to stdout.
                                      If --error-qs 0 or 2 is used, writes a single value which is used for all samples and sites.
         -printQsError [0]|1 ________ 0: Disabled, 1: Print the error probability used in quality score calculations to stdout.
                                      If --error-qs 2 is used, writes per-read quality score error probabilities to stdout.
                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.
         -printGlError [0]|1 ________ 0: Disabled, 1: Print the error probability used in genotype likelihood calculations to stdout. (requires: -GL 2)
                                      Since -GL 1 works directly with quality scores, this option is only available when -GL 2 is used.
                                      If --error-qs 2 is used, writes per-read error probabilities to stdout.
                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.
                                      If --precise-gl 1 is used, the printed values are the same as those printed by -printQsError.
         -printQScores [0]|1 ________ 0: Disabled, 1: Print the quality scores to stdout.

Output VCF/BCF tags:                  0: Do not add, 1: Add
         -addGL 0|[1] _______________ Genotype likelihoods (GL) tag
         -addGP [0]|1 _______________ Genotype probabilities (GP) tag
         -addPL [0]|1 _______________ Phred-scaled genotype likelihoods (PL) tag
         -addI16 [0]|1 ______________ I16 tag
         -addQS [0]|1 _______________ Quality score sum (QS) tag
         -addFormatDP [1]|0 _________ Per-sample read depth (FORMAT/DP) tag
         -addFormatAD [0]|1 _________ Per-sample allelic read depth (FORMAT/AD) tag
         -addFormatADF [0]|1 ________ Per-sample forward-strand allelic read depth (FORMAT/ADF) tag
         -addFormatADR [0]|1 ________ Per-sample reverse-strand allelic read depth (FORMAT/ADR) tag
         -addInfoDP [0]|1 ___________ Total read depth (INFO/DP) tag
         -addInfoAD [0]|1 ___________ Total allelic read depth (INFO/AD) tag
         -addInfoADF [0]|1 __________ Total forward-strand allelic read depth (INFO/ADF) tag
         -addInfoADR [0]|1 __________ Total reverse-strand allelic read depth (INFO/ADR) tag

</pre>
</details>

<details> <summary> Option descriptions </summary>

Option format is `-s, --long-option TYPE [X]: Description` where `-s` is the short format of the option (if any), `--long-option` is the long option, `TYPE` is the type of the argument value, and `[X]` denotes that the default argument value is `X` (if any).

Type of the argument values can be:

| Type          | Description                                            | Example                                                                |
| ------------- | ------------------------------------------------------ | ---------------------------------------------------------------------- |
| `INT`         | Integer                                                | `--long-option 3`                                                      |
| `INT+`        | Additive integer                                       | `--long-option 3` will do both `--long-option 1` and `--long-option 2` |
| `FLOAT`       | Floating point number                                  | `--long-option 3.14`                                                   |
| `STRING`      | String                                                 | `--long-option "hello"`                                                |
| `FILE`        | Filename                                               | `--long-option "file.txt"` or `--long-option /path/to/file.txt`        |
| `x\|y\|z`     | One of the listed values x, y, or z                    | `--long-option x`                                                      |
| `INT(,INT..)` | A single integer or a comma-separated list of integers | `--long-option 1,2,3` or `--long-option 1`                             |

</details>

## General Options

- `-v, --verbose INT [0]`: Verbosity level.
- `-@, --threads INT [1]`: Number of threads.
- `-s, --seed INT [time]`: Random seed for initializing the random number generator. If not defined, the seed will be set using the current system time.

  N.B. If multiple simulations are run in parallel and the seed is not defined, the random number generator may be initialized with the same seed for all of the simulations, given that it is likely that the simulations start at the same time. Therefore, it is highly recommended to set the seed manually when running multiple simulations in parallel.

## Input/Output

- `-i, --input FILE`: Input filename.
- `-o, --output STRING ['output']`: Output filename prefix. If not defined, the output filename prefix will be set to `output`. The program will append the suffix to the output filename prefix based on the given options. A simulation argument values log file (ending with `.arg`) is always generated.

  - Example: `vcfgl -i test/data/data1.vcf -o test/data1_output -e 0.1 -d 1 ` will generate `data1_output.arg`, `data1_output.bcf` files in the `test` directory.

  - Example: `vcfgl -i test/data/data1.vcf -e 0.1 -d 1 -printPileup 1 -printTruth 1 -O z -o data1` will generate `data1.arg`, `data1.pileup.gz` (`-printPileup 1`), `data1.truth.vcf.gz` (`-printTruth 1`), and `data1.vcf.gz` (`-O z`) files in the current directory.

- `-O, --output-mode [b]|u|z|v`: Output file format.
  - `b`: [Default] Compressed BCF (`.bcf`)
  - `u`: Uncompressed BCF (`.bcf`)
  - `z`: Compressed VCF (`.vcf.gz`)
  - `v`: Uncompressed VCF (`.vcf`)

## Simulation Parameters

- `-d, --depth FLOAT|'inf'`: Mean per-site read depth
- `-df, --depths-file FILE`: File containing mean per-site read depth values for each sample (one value per line)
- `-e, --error-rate FLOAT`: Base-calling error probability
- `-eq, --error-qs [0]|1|2`: Error simulation type:
  - `0`: Do not simulate errors in quality scores
  - `1`: Simulate errors in the per-site probability of wrong base calls (requires `--beta-variance`)
  - `2`: Simulate errors in the reported quality scores and genotype likelihoods (requires `--beta-variance`)
- `-bv, --beta-variance FLOAT`: Designated variance for the beta distribution
- `-GL, --gl-model 1|[2]`: Genotype likelihood model to be used in simulation
  - `1`: Genotype likelihood model with correlated errors (Li model)
  - `2`: Canonical genotype likelihood model with independent errors (McKenna model)
- `--gl1-theta FLOAT [0.83]`: Theta parameter for genotype likelihood model 1 (requires `-GL 1`)
- `--platform [0]|1`: Use platform specification
- `--precise-gl [0]|1`: Use precise error probabilities in the genotype likelihood calculation (requires `-GL 2`)
- `--i16-mapq INT [20]`: Mapping quality score for I16 tag (requires `-addI16 1`)
- `--gvcf-dps INT(,INT..)`: Minimum per-sample read depth range(s) for constructing gVCF blocks (requires `-doGVCF 1`)
  - Example: `--gvcf-dps 5,10,20` will group invariable sites into three types of gVCF blocks: `[5,10)`, `[10,20)`, and `[20,inf)`. Sites with minimum depth < 5 will be printed as regular VCF records.
- `-explode [0]|1`: Explode to sites that are not in the input file
- `--rm-invar-sites INT+`: Remove invariable sites
- `--rm-empty-sites [0]|1`: Keep or remove empty sites in the output file
- `-doUnobserved INT [1]`: Trim unobserved alleles
- `-doGVCF [0]|1`: Output in GATK GVCF format
- `-printPileup [0]|1`: Output in pileup format
- `-printTruth [0]|1`: Output the VCF file containing the true genotypes

## Output VCF/BCF Tags

Control which tags are added to the output VCF/BCF files.

- `-addGL 0|[1]`: Genotype likelihoods (GL) tag
- `-addGP [0]|1`: Genotype probabilities (GP) tag
- `-addPL [0]|1`: Phred-scaled genotype likelihoods (PL) tag
- `-addI16 [0]|1`: I16 tag
- `-addQS [0]|1`: Quality score sum (QS) tag
- `-addFormatDP [1]|0`: Per-sample read depth (FORMAT/DP) tag
- `-addFormatAD [0]|1`: Per-sample allelic read depth (FORMAT/AD) tag
- `-addFormatADF [0]|1`: Per-sample forward-strand allelic read depth (FORMAT/ADF) tag
- `-addFormatADR [0]|1`: Per-sample reverse-strand allelic read depth (FORMAT/ADR) tag
- `-addInfoDP [0]|1`: Total read depth (INFO/DP) tag
- `-addInfoAD [0]|1`: Total allelic read depth (INFO/AD) tag
- `-addInfoADF [0]|1`: Total forward-strand allelic read depth (INFO/ADF) tag
- `-addInfoADR [0]|1`: Total reverse-strand allelic read depth (INFO/ADR) tag

# Tutorials

# Contact

If you have any questions about the program, feature requests, or bug reports, you can open a new issue at [GitHub Issues](https://github.com/isinaltinkaya/vcfgl/issues). Alternatively, you can contact the main developer using the contact information at [isinaltinkaya.com](http://isinaltinkaya.com/).

# How to cite

<!-- For more details, please refer to the documentation at [](). -->

<!--

1. For each site; for each haplotype; sample number of reads from Poisson distribution with mean equal to given depth value

2. Simulating error

(2.1) Probability of error: e

(2.2) Probability of sampling an incorrect allele due to error: error rate × 1/3

(2.3) Probability of sampling correct allele: 1 − error rate

4. Normalize the log10 genotype likelihood by substracting the maximum genotype likelihood value
observed at site. -->

<!--
## Installation

To install **vcfgl**, simply download the binary executable for your platform from the [latest releases](https://github.com/your-username/vcfgl/releases) on GitHub.
### Method 0: Using release tarball

//// TODO

 -->
<!--


## Input file

Input is typically a VCF/BCF file obtained by `tskit.write_vcf`.

Input file shoud have binary haplotypes set as REF and ALT alleles. For obtaining this format from mutation simulations, you can use the binary model with `sim_mutations`:

```python
mut_ts = msprime.sim_mutations(ts, ..., model="binary")
with open("simulated.vcf","w") as vcf:
    mut_ts.write_vcf(vcf)
```

For more information about the python code please refer to <https://tskit.dev/tskit/docs/stable/python-api.html>.



## Simulate unobserved invariable sites (`-explode 1`)

You can use `-explode 1` to expand simulations to the sites not observed in the input VCF file.


Example command:
```shell
./vcfgl -i example.vcf -e 0.01 -d 5 -O v -d 3 --seed 42 -explode 1
```

Example input: `example.vcf`

Our input file contains two sites, and its header contains one chromosome named "chr22" with length 5.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=5>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
chr22	1	.	0	1	.	PASS	.	GT	0|0	0|0
chr22	3	.	0	1	.	PASS	.	GT	1|1	1|1
```


Example output: `output.vcf`

The missing 3 sites (at position 2,4, and 5) are generated by `-explode 1` by setting all individuals' genotypes to `0|0`. Thus the output file contains information for all sites.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr22,length=5>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##fileDate=Tue Sep 12 15:46:49 2023
##source=vcfgl --input example.vcf --output output --output-mode v --error-rate 0.010000 --depth 3.000000 --depths-file (null) --seed 42 -explode 1 -printBaseCounts 0 -addGP 0 -addPL 0 -addI16 0 -addQS 0
##source=vcfgl version: v0.3.0 ee51b23-dirty
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood in log10 likelihood ratio format">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Simulated per-sample read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2
chr22	1	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551	0|0:5:0,-1.49785,-12.3638,-1.49785,-12.3638,-12.3638,-1.49785,-12.3638,-12.3638,-12.3638
chr22	2	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:3:0,-0.898711,-7.41827,-0.898711,-7.41827,-7.41827,-0.898711,-7.41827,-7.41827,-7.41827	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551
chr22	3	.	A	C,G,T	.	PASS	.	GT:DP:GL	1|1:1:-2.47276,-0.29957,0,-2.47276,-0.29957,-2.47276,-2.47276,-0.29957,-2.47276,-2.47276	1|1:4:-9.89103,-1.19828,0,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103
chr22	4	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:5:0,-1.49785,-12.3638,-1.49785,-12.3638,-12.3638,-1.49785,-12.3638,-12.3638,-12.3638	0|0:2:0,-0.59914,-4.94551,-0.59914,-4.94551,-4.94551,-0.59914,-4.94551,-4.94551,-4.94551
chr22	5	.	A	C,G,T	.	PASS	.	GT:DP:GL	0|0:4:0,-1.19828,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103,-9.89103	0|0:4:0,-1.19828,-9.89103,-1.19828,-9.89103,-9.89103,-1.19828,-9.89103,-9.89103,-9.89103
```


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

___

## Notes for ANGSD/msToGlf users

The basic functionality of the program (`--error-qs 0`) is designed to simulate data equivalent to those simulated using [msToGlf](https://github.com/ANGSD/angsd/blob/master/misc/msToGlf.c).
Similar to msToGlf, vcfgl can simulate the genotype likelihoods using direct genotype likelihood model via `-GL 2`, and output pileup format via `-printPileup 1`.


```shell
msToGlf -in input.ms -out output -err 0.01 -depth 1 -pileup 1
```

is equivalent to

```shell
vcfgl -i input.vcf -o output --error-qs 0 -e 0.01 -d 1 -GL 2 -printPileup 1
```
## Errors in quality scores

You can use `--error-qs 1` to simulate errors in quality scores by simulating the site-specific errors in the probability of wrong base calls not accounted in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution-sampled error rates and the reported quality scores use the beta distribution mean (i.e. the error rate).

You can use `--error-qs 2` to simulate errors in quality scores by simulating the errors in the reported quality scores and genotype likelihoods. The error-base choosing uses the beta distribution mean (i.e. the error rate) and the reported quality scores use the beta distribution-sampled error rates.


## Simulate base qualities for a specific sequencing platform

You can use `--platform 1` to simulate base qualities for a specific sequencing platform. Currently, only NovaSeq 6000 is supported. The qualities are binned into four possible quality values: 2, 12, 23 and 37.


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
##source=vcfgl --verbose 0 --threads 1 --input test/t1.vcf --output test/t1_out --output-mode b --depth 1.000000 --error-rate 0.010000 --error-qs 0 --beta-variance -1.000000e+00 --precise-gl 1 --seed 42 --rm-invar-sites 0 --rm-empty-sites 0 -doUnobserved 1 -doGVCF 0 --platform 0 -explode 0 -addGL 1 -addGP 0 -addPL 0 -addI16 0 -addQS 0 -addFormatDP 0 -addFormatAD 0 -addFormatADF 0 -addFormatADR 0 -addInfoDP 0 -addInfoAD 0 -addInfoADF 0 -addInfoADR 0
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
```

### - Ancestral and derived alleles

By default, `vcfgl` treats the allelic state `0` in input file as `A`, and `1` as `C`. Therefore a genotype of `00` corresponds to `AA`, `01` to `AC`, `10` to `CA`, and `11` to `CC`.

Binary representation | Allele
-- | --
0 | A
1 | C

### - Sampling strands

The following options require strand information:
* `-addI16 1`
* `-addFormatADF 1`
* `-addFormatADR 1`
* `-addInfoADF 1`
* `-addInfoADR 1`

These options will trigger the random sampling of strands. The strands are sampled as forward and reverse strands from a binomial distribution with equal probability.


### - Unobserved alleles

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
