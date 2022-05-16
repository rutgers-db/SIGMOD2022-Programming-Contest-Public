#CC = /usr/bin/g++-4.8
CC = g++

OPTION = -I./ -DIN_PARALLEL -fopenmp
CFLAGS = -std=c++14 -O3 $(OPTION) -c
LFLAGS = -std=c++14 -O3 $(OPTION)

OBJS = simjoin.o SetJoin.o OvlpJoin.o CSVReader.o Table.o

all : blocking
	rm *.o

blocking : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $@

simjoin.o : simjoin.cc CSVReader.o Table.o SetJoin.o OvlpJoin.o
	$(CC) $(CFLAGS) $< -o $@ 

CSVReader.o : CSVReader.cc Table.o
	$(CC) $(CFLAGS) $< -o $@ 

Table.o : Table.cc
	$(CC) $(CFLAGS) $< -o $@ 

SetJoin.o : SetJoin.cc
	$(CC) $(CFLAGS) $< -o $@ 

OvlpJoin.o : OvlpJoin.cc
	$(CC) $(CFLAGS) $< -o $@ 

clean :
	rm blocking
	rm *.o
