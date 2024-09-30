CPP = g++ # --std=c++11
CPPFLAGS = -O2 -Wall
CPPFLAGS += $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs) 

all : run

run : calc.o Fluc.o 
	$(CPP)  -o run calc.o Fluc.o $(LDFLAGS)

.cpp.o : 
	$(CPP) $(CPPFLAGS) -c $<

clean : 
	rm *.o run
