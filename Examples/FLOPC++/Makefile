CoinDir  := $(HOME)/COIN
SymDir  := $(CoinDir)/SYMPHONY
FlopcDir := $(CoinDir)/flopc
#LPSolverDir :=$(HOME)/OSL/osllib

CXX     =  /usr/local/bin/g++

## Chose optimization level
CXXFLAGS += -g
#CXXFLAGS += -O

##############################################################################
IncDir	:= $(FlopcDir)/include 
IncDir	+= $(CoinDir)/include 
IncDir	+= $(SymDir)/include 
#IncDir += $(LPSolverDir)/include

SHLINKPREFIX = -Wl,-rpath,

CXXFLAGS += $(addprefix -I,${IncDir}) -DHAS_RANDOM -DHAS_SRANDOM

LIBDIRS := $(CoinDir)/lib
LIBDIRS += $(FlopcDir)/lib
LIBDIRS += $(SymDir)/lib
#LIBDIRS += $(LPSolverDir)/lib

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(addprefix $(SHLINKPREFIX),$(LIBDIRS))

LIBS 	+= -lOsiSym
LIBS 	+= -lsym
LIBS	+= -lflopc
LIBS 	+= -lOsi
LIBS 	+= -lOsiClp
LIBS 	+= -lClp
#LIBS 	+= -lOsiOsl
#LIBS 	+= -losl
LIBS 	+= -lCoin
LIBS 	+= -lCgl
#LIBS 	+= -lgz

LDFLAGS	+= $(LIBS)

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default

default: all

###############################################################################

all :	tap

clean:	
	@rm tap
###############################################################################
