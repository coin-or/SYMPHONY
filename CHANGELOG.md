## CHANGELOG

### Release 5.6.18
  * Switch CI to Github Actions
  * Switch to automatically generated README
  * Update documentation for coinbrew and get latex2html working again	
  * Fix reported bugs (See #171, #177, #180)
	
### Release 5.6.17
  * Update appveyor
  * Fix several reported bugs (see PRs 10-12, Issue 2)

### Release 5.6.16
  * Fix problem with appveyor configuration

### Release 5.6.15
  * Fix for configuration issue with OpenMP when building with Visual Studio
    compiler under Msys2
  * Fix for bugs in PVM version
  * Fix for bugs in computing lower bounds
  * Fix for compilation error with gcc 6
  * Other minor fixes
  * Enabling binary uploads with BinTray

### Release 5.6.14
  * Fixing small bug with re-setting environment in user applications.
  * Fixing some memory errors arising with applications when extra variables are used.
  * Fixing small bug with tracking variable indices in branching.
  * Moving code inside critical block to prevent memory access errors in shared memory parallel mode.
  * Added switches for turning hot starting on and off
  * Updates to documentation
  * Small fixes
  * Added support for Travis and Appveyor

### Release 5.6.13
  * Release to update externals and pick up bug fixes in other projects

### Release 5.6.12
  * Fixed function for determining duals and reduced costs.
  * Making it possible to build with cut validity checking enabled.
  * Fixed bug with re-using of environment for new instances.
  * Added parameter that allows saving of LP relaxations for debugging.
  * Added function for writing LP files.
  * Fixed bug that arose when we fixed a variable in strong branching (now, we keep going with LP loop when this happens).
  * Fixed long-standing bug that arose when child 0 could be pruned while generating children.

### Release 5.6.11
  * Updated externals

### Release 5.6.10
  * Fixed bug with using more the 1 process in PVM
  * Print to stderr on ctrl-c
  * Fixed double free with prep level < 1

### Release 5.6.9
  * Fixed memory leak
  * Fixed double free when nodes are retained in memory after pruning
  * Added ability to access solution pool after solve to retrieve additional suboptimal solutions.
  * Updates to documentation 

### Release 5.6.8
  * Fix to dependency linking.
  * Fix for installation with ```DESTDIR```

### Release 5.6.7
  * Fixes to distributed parallel (PVM) version.
  * Enable primal heuristics in distributed version.
  * Update externals to most recent stables.
  * Make dependency linking default.

### Release 5.6.6:

 * Disabling pre-processor for multicriteria instances.

### Release 5.6.5:

 * Added explicit dependence on libgomp, which is needed for linking with 
compilers that don't supprt OpenMP (clang on OS X)

### Release 5.6.4:

 * Fix to pkgconfig file to add flags for OpenMP.

 * Fixes for documentation.

 *  More fixes for dependency linking.

### Release 5.6.3:

 * Fixes to applications to allow some preprocessing, which is needed for
some primal heuristics to work.

 * Fixes to preprocessor settings so that the minimal amount of preprocessing
is always done.

 *  Fixes to some methods so they still work preprely even when preprocessing is
not done.

 *  Fixes to pre-processor for applications that construct the LP relaxation 
algorithmically. 

 * Fix to pkgconfig file for applications.

 * Fix for dependency linking.

 * Other small bug fixes.

### Release 5.6.2:

 * Updates and fixes to manual and documentation files.

 * Fixes for configuring with PVM.

 *  Fixes to allow dependency linking with the application library. 

 *  Bug fix for reliability branching. 

 * Bug fix for addition of column cuts

 * Updates to externals to fix bugs in dependent projects.

### Release 5.6.1:

 * Small fixes for OpenMP

### Release 5.6.0:

 * Major re-design of shared memory parallel mode for increased efficiency and stability.

 * Fixes for distributed memory parallel version (PVM)

 *  Fixes for bicriteria version

### Release 5.5.7:

 * More updates to build harness

### Release 5.5.6:

 * More updates to build harness

### Release 5.5.5:

 * More updates to build harness

### Release 5.5.4:

 * Fix memory leak
 * Delete superfluous header file
 *  More updates to build harness

### Release 5.5.3:

 * More updates to build harness

### Release 5.5.2:

 * Fix problems with Visual Studio project files
 * Update to build tools

### Release 5.5.1:

 * Fix bug that caused all user applications to crash

### Release 5.4.8:

 * Updates to MSVC++ files (applications now use property sheets and VRP app file is fixed).
 * Other fixes for build tools.

### Release 5.4.7:

 * Updates to documentation

### Release 5.5.0:

 * Improvements to preprocessing
 * Improvements to heuristics
 * Improvements to MSVC++ support
 * Bug fixes
 * Significant performance gains 

### Release 5.4.6:
 * More fixes to allow use of CPLEX as LP solver
 * Fixes to interface with GMPL

### Release 5.4.5:

 * Fixes to allow use of CPLEX as LP solver
 * Fixes to interface with GMPL

### Release 5.4.4:

 * Fixes for build system
 * Other minor fixes

### Release 5.4.3:

 * Updates to documentation.
 * Fix to allow box-constrained integer programs.
 *  Fix for GMPL integration
 *  Fix for readline versions

### Release 5.4.2:

 * Updates to MSVC++ version 10 files

### Release 5.4.1

 * Addition of MSVC++ version 9 files

### Release 5.4.0:

 * Change license to EPL.

 * Support for MSVC++ version 10 added.

 *  Support for BuildTools version 0.7 to incoorporate recent enhancements, including proper library versioning in Linux, prohibiting installation of private headers, etc.

 *  Enhancements to unit testing.

 * Updating externals to new stable versions of dependent projects.

### Release 5.3.4:

 * Fixes to the shared memory parallel version (OpenMP). It is now pretty
stable, though some minor memory conflict conditions may arise (infrequently).

 * Fixes to allow all applications to build and run properly.

 *  Updates to documentation.

### Release 5.3.3:

 * Fixes to the build system.

### Release 5.3.2:

 * Fixes to the build system.

### Release 5.3.1:

 * Fixes to the build system.

### Release 5.3.0:

 * Major changes to the build system to allow buinding against installed
binaries, provide pkg-config support, etc.

### Release 5.2.4:

 * Fixes to restore functionality of the bicriteria solution capability.

 * Fixes to examples.

### Release 5.2.3:

 * Updates to manual.

 * Added hooks to enable the use of VRPH (https://projects.coin-or.org/VRPH) within the VRP solver.

### Release 5.2.2:

 * Bug fix release.

### Release 5.2.1: 

 * Bug fix release.

### Release 5.2.0:

 * SYMPHONY has a preprocessor now.

 * Feasibility pump primal heuristic implemented.

 *  Reliability branching is now the default branching strategy.

 *  Several new statistics now part of default output.

 * Correct setting of granularity of objective function value by calculating
   GCD of coefficients.

 * Several changes in management of valid inequalities, quality checks and
   detection of duplicacy.

 * Minor changes in management of LP solver interface.

 * Several small bug-fixes and improvements.

### Release 5.1.10:

 * New dependencies.

### Release 5.1.9:

 * New dependencies.

### Release 5.1.8:

 * Introduced use of LP hot starting.

 * Improved management of cut generation.

 *  Updated externals

 *  Minor bug fixes

### Release 5.1.7:

 * Minor bug fixes

### Release 5.1.6:

 * Only a single header file (symphony.h) needs to be installed and user 
applications only need to be able to find this one header file.

 * Fixes to MSVC++ project files.

 *  Removed dependence on qsortucb routines.

### Release 5.1.5:

 * Added support for automatic download and build of Glpk (for reading of GMPL
files).

 * Minor bugs fixed and compiler warnings eliminated.

 *  Updates to MS Visual Studio files.

 *  Added short installation verification test.

### Release 5.1.4:

 * Added ability to read files in LP format.

 * Additional configuration options.

 *  Support for new classes of cutting planes.

 *  Improved algorithm control mechanism.

 * Improved output format and additional output options.

 * Improved signal handling.

 * Shared memory parallel version tested with OpenMP in Linux and Windows.

 * Added release configuration to MSVC++ build files.

 * Improved warm starting.

 * Fixes for configuration with SoPlex and Xpress.

 * Fixed configuration on PowerPC architectures.

### Release 5.1.3:

 * Support for building static executables in Unix-like environments.

 * Improved signal-catching behavior in Unix-like environments.

 *  Updated documentation.

### Release 5.1.2:

 * Update of externals.

 * Updated documentation.

### Release 5.1.1:

 * Fixes for building in the Solaris operating system.

 * Fixes for using the GNU autotools to build with the cl compiler.

 * Fixes for sym.mak file in order to allow building with MSVC++ nmake utility.

 *  Fixes for building the unit test in the MSVC++ IDE.

 * Updated documentation

### Release 5.1.0:

 * SYMPHONY now has an interactive optimizer that can be used through a
command shell. In both the sequential and parallel configurations, the user
can set parameters, load and solve instances interactively, and display
results and statistics (see below).

 * SYMPHONY now supports automatic configuration using the new COIN-OR build
system and the GNU autotools.Using autotools utilities, it is now possible to
build SYMPHONY in most operating systems and with most common compilers
compilers without user intervention.

 * Both the distributed and shared memory parallel configurations are now
fully debugged, tested, and supported. The user can now build and execute
custom SYMPHONY applications in parallel, as well as solving generic MILPs in
parallel "out of the box."

 * There are now additional options for warm starting. The user can trim the
warm starting tree before starting to resolve a problem. More specifically,
the user can decide to initiate warm starting with a predefined partition of
the final branch-and-cut tree resulting from a previous solution procedure.
This partition can include either a number of nodes created first during the
solution procedure or all of the nodes above a given level of the tree.

 * The COIN-OR repository, the current host of SYMPHONY has recently undergone 
some significant improvements of its own that have resulted in improved 
services to users. These include: 

   * SYMPHONY has a new development Web site, where users can submit trouble
     tickets, browse the source code interactively, and get up-to-date
     information on development. The address of the new site is
     https://projects.coin-or.org/SYMPHONY.

   * SYMPHONY is now hosted using subversion, a version control system with
     features vastly improved over CVS, the previous hosting software. This
     has required some reorganization and renaming of the header files.

   * SYMPHONY is now more tightly integrated with other COIN-OR projects. Due
     to improved procedures for producing stable releases, it will now be much
     easier for us to determine the exact version of SYMPHONY and all other
     COIN projects you are using when you report a bug.

   * SYMPHONY is now distributed with all COIN software needed to build a
     complete solver. Previously, other COIN softrware packages had to be
     downloaded and installed separately.

 * Two features have been deprecated and are no longer supported:

   * The native interfaces to OSL and CPLEX are now deprecated and no longer
supported. These solvers can be called through the COIN-OR OSI interface.

   * Column generation functionality has also been officially deprecated. For
     now, there are a number of other software packages that offer better
     functionality than SYMPHONY for implementing branch and price algorithms.

 * CHANGES TO THE USER INTERFACE (FROM SYMPHONY 5.0)

   * There was one minor change to the user callback API from version 5.0 to
     5.1. The user can now execute a primal heuristic in the
     user_is_feasible() callback and return the solution to SYMPHONY. The API
     for the user_is_feasible() subroutine is now
     ```C
     int user_is_feasible(void *user, double lpetol, int varnum, int *indices,
		          double *values, int *feasible, double *objval,
		          char branching, double *heur_solution)
     ```
Any feasible solution can be passed (in dense format) through the last
argument to this function.

   * Several new subroutines were added to the callable library API.

   * The name of the header file containing the SYMPHONY API has been changed
     from ```symphony_api.h``` to ```symphony.h``` (though the former has been
     retained for backword compatibility purposes).

