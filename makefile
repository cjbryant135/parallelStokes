#Makefile for the polar SOR solver
#########################################
#options
CPP = ${CXX}
OPT = -O3
CXXFLAGS = -std=c++11

polarStokesGrid.o: polarStokesGrid.hpp polarStokesGrid.cpp
	$(CPP) -c polarStokesGrid.cpp $(CXXFLAGS) $(OPT)

main.o: main.cpp 
	$(CPP) -c main.cpp $(CXXFLAGS) $(OPT)

main_exe: main.o polarStokesGrid.o
	$(CPP) -o main_exe main.o polarStokesGrid.o -lm $(CXXFLAGS) $(OPT)

clean:
	rm -r *_exe *.o *.out 
