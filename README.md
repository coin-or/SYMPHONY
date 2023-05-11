# SYMPHONY

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/SYMPHONY?sort=semver)](https://github.com/coin-or/SYMPHONY/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

SYMPHONY is an open-source generic MILP solver, callable library, and
extensible framework for implementing customized solvers for mixed-integer
linear programs (MILPs). SYMPHONY can be built in various sequential and
parallel configurations for either distributed or shared memory architectures
and can be used "out of the box" as a solver for generic mixed-integer
linear programs or customized through a wide variety of user callback
functions and control parameters. SYMPHONY has a number of advanced
capabilities stemming from the research projects discussed above, including
the ability to solve multi-objective MILPs, the ability to warm start its
solution procedure, and the ability to perform basic sensitivity analyses.
SYMPHONY has has been deployed in a variety of application areas, including
computational biology, wireless telecommunications, supply chain management,
transportation services, and air transportation. 


SYMPHONY is written in C and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The SYMPHONY development site is https://github.com/coin-or/SYMPHONY.

## CITE

Code: [![DOI](https://zenodo.org/badge/172774599.svg)](https://zenodo.org/badge/latestdoi/172774599)

Paper: http://dx.doi.org/10.1007/0-387-23529-9_5

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/SYMPHONY/actions/workflows/windows-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/SYMPHONY/actions/workflows/windows-ci.yml?query=branch%3Amaster)

[![Linux and MacOS Builds](https://github.com/coin-or/SYMPHONY/actions/workflows/linux-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/SYMPHONY/actions/workflows/linux-ci.yml?query=branch%3Amaster)

## DOWNLOAD

What follows is a quick start guide for obtaining or building
SYMPHONY on common platforms. More detailed information is
available [here](https://coin-or.github.io/user_introduction.html).

### Docker image

There is a Docker image that provides SYMPHONY, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/SYMPHONY/releases). Older binaries
are archived as part of SYMPHONY
[here](https://www.coin-or.org/download/binary/SYMPHONY).

 * *Linux* (see https://repology.org/project/coin-or-symphony/versions for a complete listing): 
   * arch:
     ```
     $ sudo pacman -S  coin-or-symphony
     ```
   * Debian/Ubuntu:
     ```
     $ sudo apt-get install  coinor-symphony coinor-libsymphony-dev
     ```
   * Fedora/Redhat/CentOS:
     ```
     $ sudo yum install  coin-or-SYMPHONY coin-or-SYMPHONY-devel
     ```
   * freebsd:
     ```
     $ sudo pkg install math/symphony
     ```
   * linuxbrew:
     ```
     $ brew install symphony
     ```
 * *Windows*: The easiest way to get SYMPHONY on Windows is to download an archive as described above.
 * *Mac OS X*: The easiest way to get SYMPHONY on Mac OS X is through [Homebrew](https://brew.sh).
     ```
     $ brew tap coin-or-tools/coinor
     $ brew install coin-or-tools/coinor/symphony
     ```

Due to license incompatibilities, pre-compiled binaries lack some 
functionality. If binaries are not available for your platform for the latest 
version and you would like to request them to be built and posted, feel free 
to let us know on the mailing list. 

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of SYMPHONY from the
 [releases](https://github.com/coin-or/SYMPHONY/releases) page,
 * Cloning this repository from [Github](https://github.com/coin-or/SYMPHONY), or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

### Dependencies

SYMPHONY has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

These quick start instructions assume you are in a bash shell. 

### Using `coinbrew`

To download and build SYMPHONY from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch SYMPHONY@master
./coinbrew build SYMPHONY
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/SYMPHONY
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## USING SYMPHONY

### Using SYMPHONY from the command line

To use SYMPHONY as a generic solver, type the executable name on the command
line, followed by one or more of the command-line switches. On the
command-line, there is one required switch---you must specify the location of
the input file by using either `-F 'filename'` (MPS file or automatic
detection with file extension) or `-L 'filename'` (LP format). If the `-D`
switch is also present, the file will be assumed to be a GMPL model file with
the data file specified after the `-D` switch. In LINUX, the following command
would solve the instance `sample.mps`

```symphony -F sample.mps```

The remaining switches are used to set SYMPHONY's native parameters on the
command line. Below is a list of these parameters. This list can also be
obtained by executng 

```symphony -h```

Note that all SYMPHONY parameters are denoted by a lowercase letter. Many
other parameters can be set with the use of a parameter file (specified with
-f). These parameters are listed in the SYMPHONY user's manual.

```
symphony [ -FL file ] [ -f parameter_file_name ]
        [ -hd ] [-a 0/1] [-b 0/1 ] [-s cands] [-l 0/1] [ -q 0/1 ] [ -r 0/1]
        [-j 0/1 ] [ -e n ] [ -i iters ] [ -t time ] [ -g gap ] [ -n nodes ]
        [ -u ub ] [ -p procs ] [ -k rule ] [ -v level ] [ -c rule ]
        [ -m max ] [ -z n ] [-o tree_out_file] [-w 0/1]


        -F model: model should be read in from file 'model'
                  (MPS format is assumed unless -D is also present)
        -L model: LP format model should be read in from file 'model'
        -D data: model is in AMPL format and data is in file 'data'
        -T dir: run test with MIPLIB3 models
        -h: help
        -f file: read parameters from parameter file 'file'
        -d: stop at first feasible solution
        -a 0/1: whether to use primal heuristics
        -b 0/1: whether to use reliability branching
        -s cands: use at most 'cands' candidates for strong branching
        -l 0/1: whether to impose a limit on strong branching time
        -q 0/1: whether or not to tighten root bounds
        -r 0/1: whether or not to do reduced cost tightening
        -j 0/1: whether or not to generate cgl cuts
        -w 0/1: whether or not to use hot starting in strong branching
        -e n: set pre-processing level to 'n'
        -i iters: allow a max of 'iters' iterations in presolve
        -t time: set wallclock time limit to 'time'
        -g gap: set gap limit to 'gap'
        -n nodes: set node limit to 'nodes'
        -u ub: use initial upper bound 'ub'
        -p procs: allow 'procs' additional threads or processors
        -k i: use node selection rule 'i'
        -v n: set verbosity to level 'n'
        -c i: use rule 'i' to compare candidates
        -m max: allow a max of 'max' cuts to enter per iteration
        -z n: set diving threshold to 'n'
        -o file: output vbc-like tree information to file 'file'
```
		
### Using the SYMPHONY interactive optimizer

To use SYMPHONY's Interactive shell, run the executable name without any
command line arguments. Then type `help` or `?` to see a list of available
commands which are as follows for this version:

```
	load      : read a problem in mps or ampl format
	solve     : solve the problem
	lpsolve   : solve the lp relaxation of the problem
	set       : set a parameter
	display   : display optimization results and stats
	reset     : restart the optimizer
	help      : show the available commands/params/options	

	quit/exit : leave the optimizer
```
	
So, if you want to load and solve an ampl/gmpl file, you will need to type
`load sample.mod sample.dat` and then `solve`. 

### Using the callable library

To use SYMPHONY as a generic callable library, compile SYMPHONY as described
above. The library that is created along with the solver itself can be linked
to using the API described in the user's manual. For examples of using the
callable library in this way, see the Examples/ subdirectory.

## DEVELOPING CUSTOM APPLICATIONS

To customize SYMPHONY by implementing the custom callback functions, simply
modify the files in the SYMPHONY/Applications/USER/ subdirectory, as described
in the user's manual and follow the compilation procedures in the file
SYMPHONY/Applications/USER/README. There are a number of sample applications
available as examples of how to do this kind of development with SYMPHONY.
These include solvers for the matching problem, the set partitioning problem
(simple and advanced versions), the vehicle routing and traveling salesman
problems, and the mixed postman problem. These applications are distributed as
separate packages and can be downloaded from http://www.branchandcut.org.
There is a white paper that guides the user through the development of the
matching solver.

## CURRENT TESTING STATUS

SYMPHONY can be used in a very large number of possible configurations and
we simply aren't able to test them all. Below is a rough idea of the testing
status of various configurations to date. If you need a certain configuration,
I would be happy to help you get it running. Please let me know.

## LP INTERFACES

**The native interfaces for OSL and CPLEX have been deprecated**
**Only LP solvers with OSI interfaces are supported**

Well tested: CPLEX, CLP

Well tested, but have some stability or other issues: GLPK

Compiled, but not well tested: SPX

## TESTED CONFIGURATIONS

### SEQUENTIAL

Sequential configurations are now automatically built and tested on Linux, OS X, and Windows using Github Actions (see above for status).

### SHARED MEMORY PARLLEL (OpenMP)

Builds and passes unit test with gcc and CLP on LINUX.

### DISTRIBUTED MEMORY PARALLEL (PVM)

Known configurations that build and pass unit test

- gcc on LINUX with PVM 3.4

### APPLICATIONS

- SYMPHONY (used as a generic MILP solver): Well tested.

- MATCH (matching): Tested, but not very extensively.

- MPP (mixed postman problem): Tested, but not very extensively.

- VRP (vehicle routing problem): Well tested.

- CNRP (capacitates network routing problem): Well tested.

- MCKP (multi criteria knapsack problem): Well tested.

- SPP (set partitioning problem): Tested, but not very extensively.

- SPP+CUTS (set partitioning problem with cutting planes): Tested, but not very 
extensively.

### CUT GENERATORS

Cut generators are supplied by the Cut Generation Library (CGL). The cut
generators that are turned on by default have been well tested. Two cut
generators that are part ofthe CGL are turned off by default because of known
issues. These are lift and project cuts and the simple rounding cuts. The
generator for Gomory cuts works well, but has somenumerical issues. We found a
few cases where the optimal solution was not found when using the Gomory cut
generator, especially in combination with CPLEX. If the solver is not
performing as it should, try turning off some of the cut generators to see if
that fixes the problem. 

## ADDITIONAL DOCUMENTATION

If you have downloaded a source distribution, LaTex source for the full
documentation is available in the SYMPHONY/Doc/ subdirectory. Quick start
guides and pointers to other on-line documentation can be found below.

The PDF version of the manual is available here:

https://coin-or.github.io/SYMPHONY/doc/SYMPHONY-5.7.1-Manual.pdf

An on-line version that is for a slightly earlier release is here:

http://coin-or.github.io/SYMPHONY/man-5.7/

## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [Discussion formum](https://github.com/coin-or/SYMPHONY/discussions)
 * [Report a bug](https://github.com/coin-or/SYMPHONY/issues/new)

### AUTHORS

SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and Laci Ladanyi
(ladanyi@us.ibm.com). Menal Guzelsoy (menal.guzelsoy@gmail.com) and Ashutosh
Mahajan (amahajan@iitb.ac.in.edu) have been instrumental in development since
version 5.0.

## ACKNOWLEDGEMENT

SYMPHONY was developed with support from
  * National Science Foundation grants CMMI-1435453, CMMI-0728011, DMI-0522796,
    DMI-0534862, DMS-9527124, and CMMI-1130914.
  * Texas ATP Grant 97-3604-010
  * Cornell University
  * Lehigh University
  * Zuse Institute Berlin
  * Research Campus Modal 'Mathematical Optimization and Data Analysis 
Laboratories' funded by the German Federal Ministry of Education and Research
(BMBF Grant 05M14ZAM) and by the DFG SFB/Transregio 154

