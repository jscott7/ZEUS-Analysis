#
# Flags for cpp compile fortran
#
# DEBUGFLAGS := -DDEBUG
DEBUGFLAGS := #
# EXTRAFLAGS := -DMC -DKSHORTNT -DDSUBESSNT
# EXTRAFLAGS := -DKSHORTNT -DDSUBESSNT
EXTRAFLAGS := #
#
# Flags for compliation
#
CPPFLAGS := -cpp $(DEBUGFLAGS) $(EXTRAFLAGS)
FFLAGS   = -C
CFLAGS   = -Olimit 2000 -c -O
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

default:: ana
.PHONY: clean

#main.o: $(CSRCS)
#	$(CC) $(CFLAGS) $^ -o $(OBJDIR)/main.o

ana: 
	gfortran $(CPPFLAGS) $(FFLAGS) $(FSRCS) $(FPPRCS) $(INCRCS) \
		-o $(TARGETDIR)/ana 	\
        -L /usr/lib/i386-linux-gnu \
		-lmathlib -lpacklib -lkernlib -lpdflib804


clean:
	$(RM) $(TARGETDIR)/ana $(OBJDIR)/* $(JUNK)

