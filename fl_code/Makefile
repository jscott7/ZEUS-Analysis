# Make ZEUS Analysis code
# Run: make analysis
#
# Flags for compile fortran
#
# -fcheck=all : Add runtime checks, e.g for array bounds
# -g : Add extra debug information
# DEBUGFLAGS := -fcheck=all -g
EXTRAFLAGS := #
#
# Flags for compliation
#
# -cpp : Enable preprocessing
CPPFLAGS := -cpp $(DEBUGFLAGS) $(EXTRAFLAGS)
FFLAGS   = -O3 -C
# -c : Compile to obj files -O optimimations -O3 higher level of optimisation
# There's no real tangible difference for a single '96 batch between -O and -O3
#CFLAGS   = -Olimit 2000 -c -O3
#
# Define directories
#
TARGETDIR  := ./exe
USERSRCS   := ./src
INCLUDEDIR := ./inc
OBJDIR     := ./obj
CPPFLAGS   += -I $(INCLUDEDIR)
VPATH      := $(USERSRCS):$(TARGETDIR)

FSRCS  := $(foreach dir,$(USERSRCS),$(wildcard $(dir)/*.f))
FPPRCS  := $(foreach dir,$(USERSRCS),$(wildcard $(dir)/*.fpp))
INCRCS  := $(foreach dir,$(USERSRCS),$(wildcard $(dir)/*.inc))
CSRCS  := $(foreach dir,$(USERSRCS),$(wildcard $(dir)/*.c))
JUNK   := $(foreach dir,$(USERSRCS),$(wildcard $(dir)/*~))

default:: analysis

# PHONY ensures that make clean will run the recipe regardless of whether there is a file
# named clean
.PHONY: clean

analysis: 
	gfortran $(CPPFLAGS) $(FFLAGS) $(FSRCS) $(FPPRCS) $(INCRCS) \
		-o $(TARGETDIR)/analysis.exe 	\
        -L /usr/lib/i386-linux-gnu \
		-lmathlib -lpacklib -lkernlib -lpdflib804


clean:
	$(RM) $(TARGETDIR)/ana $(OBJDIR)/* $(JUNK)


