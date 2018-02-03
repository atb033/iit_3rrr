VPATH =  ../SWZ 

CC=g++
LD=g++ -fopenmp
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g -std=c++11 -O3 

OBJS = main.o  inertiamodel_3rrr.o matrices.o torque.o torquetest.o dynindices.o scanner.o swz_3rrr.o swz_len.o
DEPS = inertiamodel_3rrr.h scanner.h swz_len.h swz_3rrr.h

MAIN = out

$(MAIN):$(OBJS)
	$(LD) -pg $(OBJS) $(CFLAGS) -o $(MAIN) 

%.o: %.cpp $(DEPS)
	$(LD) $(CFLAGS) -c $< -pg

clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res
