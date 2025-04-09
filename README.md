<a name="readme-top"></a>

<h3 align="center">vcfgl</h3>

<p align="center">
<a href="https://github.com/isinaltinkaya/vcfgl/releases/latest"><img src="https://img.shields.io/github/v/release/isinaltinkaya/vcfgl?labelColor=white&color=blue"/></a>
  <a href="https://github.com/isinaltinkaya/vcfgl/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-GPLv3.0-purple.svg?labelColor=white"/></a>
<img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/isinaltinkaya/vcfgl/test.yml?labelColor=white">
  <a href="https://bio.tools/vcfgl"><img src="https://img.shields.io/badge/bio.tools-vcfgl-orange?labelColor=white"/></a>
<a href="https://doi.org/10.1101/2024.04.09.586324"><img src="https://img.shields.io/badge/bioRxiv-10.1101%2F2024.04.09.586324-black?labelColor=white&color=darkred"/></a>
<!-- <a href="https://github.com/isinaltinkaya/vcfgl/actions/workflows/test.yml"><img src="https://github.com/isinaltinkaya/vcfgl/actions/workflows/test.yml/badge.svg" /></a>  -->
</p>


  <p align="center">
    Genotype likelihood simulator for VCF/BCF files
    <br />
    <!-- <a href="https://github.com/isinaltinkaya/vcfgl"><strong>Quickstart»</strong></a> -->
    <br />
    <br />
    <a href="https://github.com/isinaltinkaya/vcfgl">Installation</a>
    <b>·</b>
    <a href="https://github.com/isinaltinkaya/vcfgl/doc/tutorial.MD">Tutorial</a>
    <b>·</b>
    <a href="https://github.com/isinaltinkaya/vcfgl/issues/new?assignees=&labels=bug&projects=&template=report-a-bug.md&title=%5BBUG%5D">Report Bug</a>
    <b>·</b>
    <a href="https://github.com/isinaltinkaya/vcfgl/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.md&title=%5BFR%5D">Request Feature</a>
    <b>·</b>
    <a href="#how-to-cite">Cite</a>
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
    <!-- <ol>
      <li><a href="#general-options">General Options</a></li>
      <li><a href="#inputoutput">Input/Output</a></li>
      <li><a href="#simulation-parameters">Simulation Parameters</a></li>
      <li><a href="#output-vcfbcf-tags">Output VCF/BCF Tags</a></li>
    </ol> -->
    </li>
    <li> <a href="#tutorials">Tutorials</a></li>
    <li><a href="#quickstart-for-mstoglf-users">Quickstart for msToGlf users</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#how-to-cite">How to cite</a></li>
  </ol>
</details>

# Overview

**vcfgl** is a lightweight command-line program for simulating VCF/BCF and gVCF files. It allows you to simulate sequencing data with various parameters, such as read depth, base-calling error rates, quality score errors, and genotype likelihood models.

# Installation

You can install **vcfgl** using one of the following methods:

### &rarr; Method 1: Using release tarball (recommended)

You can download the latest release tarball from the [GitHub releases page](https://github.com/isinaltinkaya/vcfgl/releases/latest).

### &rarr; Method 2: Using HTSlib submodule

This method uses the htslib submodule included in the repository.

```shell
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl;
make;
```

### &rarr; Method 3: Using systemwide HTSlib installation

This method assumes you have a systemwide htslib installation.

```shell
git clone https://github.com/isinaltinkaya/vcfgl.git;
cd vcfgl;
make HTSSRC="systemwide";
```

### &rarr; Method 4: Using specified HTSlib path

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
<pre>

<pre>

Usage: vcfgl -i <input> [options]

    -h, --help _____________________  Print this help message and exit
    -v, --version __________________  Print version and build information and exit

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
    -o, --output STRING ['output'] __ Output filename prefix
    -O, --output-mode [b]|u|z|v _____ b: Compressed BCF (.bcf), u: uncompressed BCF (.bcf), z: compressed VCF (.vcf.gz), v: uncompressed VCF (.vcf)

Simulation parameters:
    -d, --depth FLOAT|'inf' _________ Mean per-site read depth
                                      ('inf') Simulate true values (requires: -addFormatDP 0 -addInfoDP 0)
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
         --qs-bins FILE _____________ File containing the quality score binning to be used in the simulation
         --precise-gl [0]|1 _________ 0: Use the discrete phred-scaled error probabilities in the genotype likelihood calculation
                                      1: Use precise error probabilities in the genotype likelihood calculation (requires: -GL 2)
         --i16-mapq INT [20] ________ Mapping quality score for I16 tag (requires: -addI16 1)
         --gvcf-dps INT(,INT..) _____ Minimum per-sample read depth range(s) for constructing gVCF blocks (requires: -doGVCF 1)
                                      Example: `--gvcf-dps 5,10,20` will group invariable sites into three types of gVCF blocks: [5,10), [10,20) and [20,inf)
                                      Sites with minimum depth < 5 will be printed as regular VCF records.
         --adjust-qs INT+ [0] _______ 0: Do not adjust quality scores
                                      1: Use adjusted quality scores in genotype likelihoods (requires: --precise-gl 0)
                                      2: Use adjusted quality scores in calculating the quality score sum (QS) tag (requires: -addQS 1)
                                      4: Use adjusted quality scores in pileup output (requires: --printPileup 1)
                                      8: Use adjusted quality scores in --printQScores output (requires: --printQScores 1)
                                      16: Use adjusted quality scores in --printGlError output (requires: --printGlError 1)
         --adjust-by FLOAT [0.499] __ Adjustment value for quality scores (requires: --adjust-qs > 0)
         -explode [0]|1 _____________ 1: Explode to sites that are not in input file.
                                      Useful for simulating invariable sites when the input file only contains variable sites.
                                      Sets all genotypes in exploded sites to homozygous reference.
         --rm-invar-sites INT+ [0]___ 0: Do not remove invariable sites
                                      1: Remove sites where all individuals' true genotypes in the input file are homozygous reference
                                      2: Remove sites where all individuals' true genotypes in the input file are homozygous alternative
                                      4: Remove sites where the all simulated reads among all individuals are the same base
                                      Example: '--rm-invar-sites 3' (1+2) will do both 1 and 2 (i.e. remove all homozygous sites)
         --rm-empty-sites [0]|1 _____ 0: Do not remove empty sites
                                      1: Remove empty sites (i.e. sites where no reads were simulated)

Output options:
         -doUnobserved INT [1] ______ 0: Trim unobserved alleles. Only alleles that are observed will be listed in REF and ALT fields
                                      1: Use '<*>' notation to represent unobserved alleles
                                      2: Use '<NON_REF>' notation to represent unobserved alleles (a.k.a. GATK notation)
                                      3: Explode unobserved bases from {A,C,G,T} list
                                      4: Use '<*>' notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list
                                      5: Use '<NON_REF>' notation to represent unobserved alleles and explode unobserved bases from {A,C,G,T} list
         -doGVCF [0]|1 ______________ 0: Disabled
                                      1: Output in gVCF format (requires: --rm-invar-sites 0, -doUnobserved 2, -addPL 1 and --gvcf-dps INT)

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

Additional output formats:
         -printPileup [0]|1 _________ 0: Disabled, 1: Output in pileup format (<output_prefix>.pileup.gz)
         -printTruth [0]|1 __________ 0: Disabled, 1: Output the VCF file containing the true genotypes (<output_prefix>.truth.vcf)
         -printBasePickError [0]|1 __ 0: Disabled, 1: Print the base picking error probability to stdout.
                                      If --error-qs 1 is used, writes per-read base picking error probabilities to stdout.
                                      If --error-qs 0 or 2 is used, writes a single value which is used for all samples and sites.
                                      The columns are: type, sample_id, contig, site, read_index, base_pick_error_prob
         -printQsError [0]|1 ________ 0: Disabled, 1: Print the error probability used in quality score calculations to stdout.
                                      If --error-qs 2 is used, writes per-read quality score error probabilities to stdout.
                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.
                                      The columns are: type, sample_id, contig, site, read_index, error_prob
         -printGlError [0]|1 ________ 0: Disabled, 1: Print the error probability used in genotype likelihood calculations to stdout. (requires: -GL 2)
                                      Since -GL 1 works directly with quality scores, this option is only available when -GL 2 is used.
                                      If --error-qs 2 is used, writes per-read error probabilities to stdout.
                                      If --error-qs 0 or 1 is used, writes a single value which is used for all samples and sites.
                                      If --precise-gl 1 is used, the printed values are the same as those printed by -printQsError.
                                      The columns are: type, sample_id, contig, site, read_index, error_prob
         -printQScores [0]|1 ________ 0: Disabled, 1: Print the quality scores to stdout.
                                      The columns are: type, sample_id, contig, site, read_index, qScore

</pre>
</pre>
</details>




<details closed> <summary> <b> Option descriptions </b> </summary>

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

<details closed> <summary> <b> General options </b> </summary>

- `-v, --verbose INT [0]`: Verbosity level.
- `-@, --threads INT [1]`: Number of threads.
- `-s, --seed INT [time]`: Random seed for initializing the random number generator. If not defined, the seed will be set using the current system time.

  N.B. If multiple simulations are run in parallel and the seed is not defined, the random number generator may be initialized with the same seed for all of the simulations, given that it is likely that the simulations start at the same time. Therefore, it is highly recommended to set the seed manually when running multiple simulations in parallel.

</details>

<details closed> <summary> <b> Input/Output </b> </summary>

- `-i, --input FILE`: Input filename.
- `-o, --output STRING ['output']`: Output filename prefix. If not defined, the output filename prefix will be set to `output`. The program will append the suffix to the output filename prefix based on the given options. A simulation argument values log file (ending with `.arg`) is always generated.

  - Example: `vcfgl -i test/data/data1.vcf -o test/data1_output -e 0.1 -d 1 ` will generate `data1_output.arg`, `data1_output.bcf` files in the `test` directory.

  - Example: `vcfgl -i test/data/data1.vcf -e 0.1 -d 1 -printPileup 1 -printTruth 1 -O z -o data1` will generate `data1.arg`, `data1.pileup.gz` (`-printPileup 1`), `data1.truth.vcf.gz` (`-printTruth 1`), and `data1.vcf.gz` (`-O z`) files in the current directory.

- `-O, --output-mode [b]|u|z|v`: Output file format.
  - `b`: [Default] Compressed BCF (`.bcf`)
  - `u`: Uncompressed BCF (`.bcf`)
  - `z`: Compressed VCF (`.vcf.gz`)
  - `v`: Uncompressed VCF (`.vcf`)

</details>


<details closed> <summary> <b> Simulation parameters </b> </summary>


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
- `--qs-bins FILE`: File containing the quality score binning to be used in the simulation
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
- `-printBasePickError [0]|1`: Print the base picking error probability to stdout
- `-printQsError [0]|1`: Print the error probability used in quality score calculations to stdout. The columns are: type, sample_id, contig, site, read_index, error_prob
- `-printGlError [0]|1`: Print the error probability used in genotype likelihood calculations to stdout. The columns are: type, sample_id, contig, site, read_index, error_prob
- `-printQScores [0]|1`: Print the quality scores to stdout. The columns are: type, sample_id, contig, site, read_index, qScore

</details>

<details closed> <summary> <b> Output VCF/BCF tags </b> </summary>


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

</details>

# Tutorials

  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/install.MD">Installation</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/depth.MD">Simulating read depths</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/error_qs.MD">Simulating quality score errors</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/simulate_unobserved.MD">Simulating unobserved sites</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/qs_binning.MD">Simulate quality score binning</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/with_msprime.MD">Using vcfgl with msprime</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/with_stdpopsim.MD">Using vcfgl with stdpopsim</a></li>
  <li><a href="https://github.com/isinaltinkaya/vcfgl/blob/main/doc/with_SLiM.MD">Using vcfgl with SLiM</a></li>

# Quickstart for msToGlf users

The program's basic functionality (i.e., simulation without errors in quality scores, `--error-qs 0`) is designed to simulate data equivalent to those simulated using [msToGlf](https://github.com/ANGSD/angsd/blob/master/misc/msToGlf.c).
Similar to msToGlf, vcfgl can simulate the genotype likelihoods using direct genotype likelihood model via `-GL 2`, and output pileup format via `-printPileup 1`.


```
msToGlf -in input.ms -out output -err 0.01 -depth 1 -pileup 1
```

is equivalent to

```shell
vcfgl -i input.vcf -o output --error-qs 0 -e 0.01 -d 1 -GL 2 -printPileup 1
```

# Contact

If you have any questions about the program, feature requests, or bug reports, you can open a new issue at [GitHub Issues](https://github.com/isinaltinkaya/vcfgl/issues). Alternatively, you can contact the main developer using the contact information at [isinaltinkaya.com](http://isinaltinkaya.com/).

# How to cite

vcfgl is published in the Bioinformatics journal, and freely available through the following link: [Publication](http://dx.doi.org/10.1093/bioinformatics/btaf098)

You can use the BibTex entry below for referencing this program in your work:

```BibTex
@ARTICLE{Altinkaya2025-ry,
  title     = "vcfgl: a flexible genotype likelihood simulator for {VCF}/{BCF}
               files",
  author    = "Altinkaya, Isin and Nielsen, Rasmus and Korneliussen, Thorfinn
               Sand",
  journal   = "Bioinformatics",
  publisher = "Oxford University Press (OUP)",
  volume    =  41,
  number    =  4,
  pages     = "btaf098",
  month     =  mar,
  year      =  2025,
  language  = "en"
}
```



