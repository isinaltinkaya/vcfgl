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


