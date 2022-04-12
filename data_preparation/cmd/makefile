#____________________________________________________________________
#
#   Makefile for EAZE
#
#	These commands are used for generating the EAZE program.
#
#  This is the highest level file:  "makefile".
#  This file first optionally includes the "Makefile.preamble" file,
#  if it exists; this file should contain
#  the local definitions of the flags to pass to the C-PreProcessor, the 
#  compiler and the linker, 
#  and it should specify the libraries to include in linking, and any
#  local redefinitions of directories and things.
#  This may contain machine-dependent definitions, or it
#  may be a link to a machine-specific file, if these things are so
#  different.
#
#  This file then includes the global file "basic.make", which contains the
#  basic setup for the standard directory hierarchy for a project.
#  It also adds the standard directories to VPATH for finding the source files
#  and generated source files.  Finally, it defines the standard Targets
#  for cleaning up the project.
#
#  Next, this makefile includes the the "machine_dep" file of machine 
#  dependencies, taken from the Make_mgr directory; 
#  there should need to be only one of these per machine.
#  The machine specification is via the ZEUS-specific
#  environment variable ZOSTYPE, 
#  which is appended to the filename as the extension:
#  "machine_dep.$(ZOSTYPE)".
#  Similarly, the "default_rules" file of default rules is included from
#  the Make_mgr directory:  "default_rules.$(ZOSTYPE)".
#
#  The Operating System type and the machine Architecture type are necessary
#  to maintain projects on disk space shared between different machines.
#  The two environment variables ZOSTYPE and ZARCH_TYPE ZEUS-specific,
#  and they are defined in the c-shell file "setup_ZSOFT.csh", 
#  which can be sourced from the user's ".cshrc" file.
#  A similar Bourne-shell file can be used.
#
#  This file next includes  "$(ZMAKEDIR)/ZEUS/setup_zeus.make", which
#  prepares the ZEUS Make environment.  This includes the definitions of
#  a consistent set of version numbers for all the ZEUS projects and
#  a set of Macros for each project pointing to the appropriate directories
#  of the appropriate version.  These are necessary for specifying the
#  search path for include files to the C-PreProcessor.  The special Macro
#  WINDOWSINCS is defined as a list of "-I" specifications to access all
#  the ZEUS Adamo-Windows.
#
#  Next, this top-level makefile includes the makefile dependencies
#  description of the specific program:  "project_dep.make".
#  This file must contain the dependencies of the project.
#  It must define a list of subdirectories of the source files, LIB_SUBSRCS,
#  and a similar list of subdirectories of the include files, LIB_SUBINCS.
#  The source files are searched for automatically from within these
#  subdirectories.
#  
#  For the official libraries, this file adds the ZEUS MAKE Rules,
#  "$(ZMAKEDIR)/library_rules.$(ZOSTYPE)" for the given machine.
#  This is necessary since these have rules which use a Macro, FULLLIBPATH,
#  which must be predefined in the "project_dep.make" file.
#
#  For flexibility to handle special cases, this file lastly includes an
#  optional file, "Makefile.postamble".
#  This file should contain any special redefinitions of Macros or rules
#  which are necessary to get around the Standard Conventions.
#  This need not exist.
#
#  The default target of this makefile is, by convention, "help".
#  The help target should be included in the makefile dependencies file,
#  "project_dep.make",
#  and it should describe the various targets that are available.
#  
#  Then, you should be able to execute "make all".
#  
#	Dr. John J. Ryan
#	21-April-1996
#____________________________________________________________________________
#
default ::  help
#
#  For generic Unix systems, which have defined the 
#  two environment variables OSTYPE and ZARCH_TYPE,
#  we can use these lines without further machine dependence:
#
#____________________________________________________________________
#
#     Optional Preamble setups (for CPP flags, etc.)
#____________________________________________________________________
#
-include  Makefile.preamble
#____________________________________________________________________
#
#     GNU Make basic setup for directory hierarchy
#____________________________________________________________________
#
include  $(ZMAKEDIR)/basic.make
#____________________________________________________________________
#
#     GNU Make setup for Handling Machine Dependence
#____________________________________________________________________
#
include  $(ZMAKEDIR)/machine_dep.make
#____________________________________________________________________
#
#     Prepare the ZEUS environment
#____________________________________________________________________
#
include  $(ZMAKEDIR)/ZEUS/setup_zeus.make
#____________________________________________________________________
#
#     Program Dependencies
#____________________________________________________________________
#
include  eaze_job.make
#____________________________________________________________________
#
#     Add the ZEUS MAKE Rules
#____________________________________________________________________
#
include  $(ZMAKEDIR)/default_rules.make
#____________________________________________________________________
#
#     Modify things locally
#____________________________________________________________________
#
-include  Makefile.postamble
#____________________________________________________________________
#
