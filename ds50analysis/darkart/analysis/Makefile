#TARGET = analysis_jan2014note
#TARGET = ds50s2studies
TARGET = merge
#DEP=$(TARGET:.out=.d)
#SRCS = analysis_jan2014note.C
#SRCS = ds50s2studies.C
SRCS = merge.C
DICTNAME = DarkArtDict
#DARKARTDIR = ../darkart/Products
#DARKARTHEADERS = $(DARKARTDIR)/EventInfo.hh $(DARKARTDIR)/Channel.hh $(DARKARTDIR)/Baseline.hh $(DARKARTDIR)/Pmt.hh $(DARKARTDIR)/WaveformInfo.hh $(DARKARTDIR)/Pulse.hh $(DARKARTDIR)/PulseParam.hh $(DARKARTDIR)/Roi.hh $(DARKARTDIR)/Spe.hh
OBJS = $(addsuffix .o, $(notdir $(basename $(SRCS))))
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs) -lProof #-lMemStat
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)
ROOTCINT = $(ROOTSYS)/bin/rootcint
CXXFLAGS   = $(ROOTCFLAGS) -Wall -fPIC -g #-pthread 
CXXLIBS    = $(ROOTLIBS)
CC = gcc 


all: $(DICTNAME).so $(TARGET) 
#all: $(TARGET)

#%: %.o
#	$(CC) $(CXXFLAGS) $(CXXLIBS) $< -o $@
%: %.o
	$(CC) $(CXXFLAGS) $(CXXLIBS) $(DICTNAME).so $< -o $@

$(OBJS): $(SRCS)
	$(CC) $(CXXFLAGS) -c $^


$(DICTNAME).cxx: 
	$(ROOTCINT) -f $@ -c -p $(SRCS) LinkDef.h

$(DICTNAME).so: $(DICTNAME).cxx
	$(CC) -shared -o $@ $(ROOTLDFLAGS) $(CXXFLAGS) $(CXXLIBS) -I$(ROOTSYS)/include $^
	#$(CC) -shared -o $@ $(ROOTLDFLAGS) $(CXXFLAGS) $(CXXLIBS) -I/ds50/app/products/gcc/v4_8_1/Linux64bit+2.6-2.12/include/c++/4.8.1/tr1 $^
clean:
	rm  -r $(OBJS) $(DEP) $(SOBJS)  $(DICTNAME).so $(DICTNAME).cxx $(DICTNAME).h

