# Makefile for compiling NSGA-II source code
VPATH =  ../SWZ ../Dynamics

CC=g++
LD=g++ -fopenmp
RM=rm -f
CFLAGS=-Wall -ansi -pedantic  -std=c++11 -O3 

OBJS =  nsga2r.o matrices.o inertiamodel_3rrr.o torque.o dynindices.o\
		allocate.o auxiliary.o crossover.o crowddist.o decode.o display.o dominance.o eval.o \
		fillnds.o initialize.o list.o merge.o mutation.o problemdef.o rand.o rank.o\
		report.o sort.o tourselect.o\
		scanner.o swz_len.o swz_3rrr.o 
DEPS  = global.h rand.h inertiamodel_3rrr.h swz_3rrr.h swz_len.h scanner.h 
MAIN= optimize

$(MAIN):$(OBJS)
	$(LD) $(OBJS) $(CFLAGS) -o $(MAIN) 

%.o: %.cpp $(DEPS)
	$(LD) $(CFLAGS) -c $< 

clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res
	$(RM)  *.o




