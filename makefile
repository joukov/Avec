######################################################################
# Copyright 2000 - 2006 Bruno Wittmer
# This file is part of the Avec package.
#
# Avec is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Avec is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Avec; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#######################################################################
#---------------------------------------------------

OPTIFLAGS= -O3 -fomit-frame-pointer -ffast-math -funroll-all-loops -fstrength-reduce -fthread-jumps -frerun-cse-after-loop -frerun-loop-opt -fgcse -fexpensive-optimizations

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

CXXFLAGS= -c -Wall $(ROOTCFLAGS) -D__AVECROOT__ -fPIC

LDFLAGS       = 

SOFLAGS       = -shared

LD            = g++

LIBS          = $(ROOTLIBS)

HDRS          = Avec.h Avec2D.h

SRCS          = Avec.cc Avec2D.cc AvecDict.C

OBJS          = Avec.o Avec2D.o AvecDict.o

SHARED_OBJS   = Avec.so

PROGRAM       = tstavec

all:            $(SHARED_OBJS)

$(PROGRAM):     tstavec.o $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) tstavec.o $(OBJS) $(LIBS) -o $(PROGRAM)
		@echo "done"

clean:;		rm -f $(OBJS) $(PROGRAM) AvecDict.* *.so core

###
Avec.o: Avec.cc Avec.h
	g++ $(CXXFLAGS) $(OPTIFLAGS) Avec.cc

Avec2D.o: Avec2D.cc Avec2D.h
	g++ $(CXXFLAGS) $(OPTIFLAGS) Avec2D.cc

AvecDict.C: Avec.h LinkDef.h
	@echo "Generating dictionary ..."
	rootcint -f AvecDict.C -c -p Avec.h Avec2D.h LinkDef.h

Avec.so:	$(OBJS)
		$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

doc: Avec.h Avec2D.h Doxyfile
	doxygen
#---------------------------------------------------
