##############################################################################################
# export coast dir
#--------------------------------------------------------------------------------------------

export COAST_DIR=/home/marcos/Documents/corsika/corsika-77550

#--------------------------------------------------------------------------------------------

LDFLAGS   = ${CXXFLAGS} -fPIC -ggdb3
LDFLAGS  += $(shell root-config --libs)
LDFLAGS  += -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown
LDFLAGS  += -lCorsikaFileIO
LDFLAGS  += -lCorsikaIntern

CPPFLAGS  = ${CXXFLAGS} -c -fPIC -ggdb3 
CPPFLAGS += -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES  = $(wildcard *.cc)
OBJECTS   = $(patsubst %.cc, %.o, ${SCRFILES})

EXE       = footprint
all: ${EXE}

#--------------------------------------------------------------------------------------------

footprint:      footprint.o
		${CXX} $^ ${LDFLAGS} -o $@

clean:
	@rm -f *.o *.so *.png
