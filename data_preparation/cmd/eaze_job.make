#__________________________________________________________________________
#
#   Makefile for EAZE library and Program Generation
#__________________________________________________________________________
#
.PHONY : all
all : exe
#____________________________________________________________________
#
#     EAZE executable
#____________________________________________________________________
#
PROGRAM_MAIN  = eaze
PROGRAM_NAME  = jps_ntuple
PROGRAM_TYPE  = eaze
#____________________________________________________________________
#
#     Program Dependencies
#____________________________________________________________________
#
#   Define the subdirectories of the Include tree
#   ---------------------------------------------
LIB_SUBINCS  := #
USERINCS := analysis v4.6

#   Define the subdirectories of the Source tree
#   --------------------------------------------
LIB_SUBSRCS := #

USERSRCS := v4.6 programs analysis
#____________________________________________________________________
#
#    Add required include directories 
#____________________________________________________________________
#
FPPFLAGS +=   $(WINDOWSINCS)
FPPFLAGS +=   $(ZEUSINCS)
#FPPFLAGS +=   $(ZGANAINCS)
#FPPFLAGS +=   $(PHANTOMINCS)
#FPPFLAGS += -I$(CCINCDIR)
#FPPFLAGS +=   $(GEANTINCS)
#FPPFLAGS += -I$(O1INCDIR)
#FPPFLAGS += -I$(ZRINCDIR)/zrcmn
#FPPFLAGS += -I$(VCINCDIR)/vckeep
#
#_____________________________________________________________________________
#
#    Use the standard GNU Make tools for EAZE Linking
#_____________________________________________________________________________
#
include  $(ZMAKEDIR)/general_build.make
include  $(ZMAKEDIR)/eaze_run.make
#
WINDOWSVERS := 1998a
#____________________________________________________________________
#

testlib : 
	$(SAY) "The project library is:  $(FULLLIBPATH)."
	$(SAY) "  "
	$(SAY) "The Fortran Source files are:"
	$(SAY) "$(SOURCESF)."
	$(SAY) "  "
	$(SAY) "The C-type Source files are:"
	$(SAY) "$(SOURCESC)."
	$(SAY) "  "
	$(SAY) "The compiled object modules are:"
	$(SAY) "$(OBJ_FILES)."
	$(SAY) "  "

testexe : 
	$(SAY) "The program name = $(PROGRAM_NAME)."
	$(SAY) "  "
	$(SAY) "The Fortran Source files are:"
	$(SAY) "$(USER_SOURCESF)."
	$(SAY) "  "
	$(SAY) "The C-type Source files are:"
	$(SAY) "$(USER_SOURCESC)."
	$(SAY) "  "
	$(SAY) "The compiled object modules are:"
	$(SAY) "$(USER_OBJ_FILES)."
	$(SAY) "  "

help : 
	$(SAY) "This procedure makes an example EAZE job."
	$(SAY) "It books and fills an example Ntuple."
	$(SAY) "It also extracts a small number of events to an output file."
	$(SAY) "  "
	$(MAKE) help_eaze
#____________________________________________________________________
#

