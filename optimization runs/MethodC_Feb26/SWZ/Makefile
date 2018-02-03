
CC=g++
LD=g++ -fopenmp
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g -std=c++11 -O3 

OBJS = main.o swz_len.o scanner.o swz_3rrr.o 
DEPS = swz_3rrr.h swz_len.h scanner.h

MAIN = out

$(MAIN):$(OBJS)
	$(LD) -pg $(OBJS) $(CFLAGS) -o $(MAIN) 

%.o: %.cpp $(DEPS)
	$(LD) $(CFLAGS) -c $< -pg

clean:
	$(RM) *.o stats timings roots testsuite *.bak *~ core *.res
	$(RM)  *.o
