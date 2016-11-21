# ilp fvs -- ilp formulation for fvs
# Copyright (C) 2014  Falk Hüffner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# configuration

HAVE_CPLEX=1
HAVE_GUROBI=1

OPTDIR=/opt
#OPTDIR=/afs/akt.tu-berlin.de/user/hueffner/opt

CPLEXVERSION=1263
CPLEXARCH=x86-64_linux
CPLEXDIR=$(OPTDIR)/cplex$(CPLEXVERSION)/cplex

GUROBIVERSION=65
GUROBIMINORVERSION=0
GUROBIARCH=linux64
GUROBIDIR=$(OPTDIR)/gurobi$(GUROBIVERSION)$(GUROBIMINORVERSION)/$(GUROBIARCH)

# end configuration

EXECS	  = fvs
LP_OBJ	  = LP/LP.o

CXX	  = g++
CXXFLAGS  = -std=c++11 -g -Ofast -march=native -W -Wall -Werror # -D_GLIBCXX_DEBUG
# CXXFLAGS += -DLP_SINGLE_THREAD
DEFINES	  =
INCLUDES  =
LIBS	  =

ifdef HAVE_CPLEX
DEFINES	  += -DHAVE_CPLEX
LP_OBJ	  += LP/LPCplex.o
INCLUDES  += -I$(CPLEXDIR)/include
LIBS	  += -L$(CPLEXDIR)/bin/$(CPLEXARCH)	\
	     -lcplex$(CPLEXVERSION)		\
	     -ldl
endif
ifdef HAVE_GUROBI
DEFINES	  += -DHAVE_GUROBI
LP_OBJ	  += LP/LPGurobi.o
INCLUDES  += -I$(GUROBIDIR)/include
LIBS      += -L$(GUROBIDIR)/lib			\
             -lgurobi$(GUROBIVERSION)
endif

all: $(EXECS)

fvs: fvs.o Util.o Graph.o $(LP_OBJ)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(LIBS) $^ -o $@

%.o: %.cc
	@mkdir -p .deps/LP
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDES) -c -MD -o $@ $<
	@cp $*.d .deps/$*.P;					    \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .deps/$*.P;    \
		rm -f $*.d

clean:
	rm -f *.o LP/*.o $(EXECS) core gmon.out

-include $(wildcard .deps/*.P)
