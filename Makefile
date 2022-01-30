all: Output

Output: main.o FractionalStepSim.o fileReadingFunctions.o fractionalStepGrid.o  general_computation_functions.o grid.o testing_functions.o multigrid.o FracStepMultigrid.o 
	g++ -Ofast -I/home/anand/ main.o FractionalStepSim.o fileReadingFunctions.o fractionalStepGrid.o  general_computation_functions.o grid.o testing_functions.o multigrid.o FracStepMultigrid.o -o Output

main.o: main.cpp
	g++ -Ofast -I/home/anand/ -c main.cpp 

FractionalStepSim.o: FractionalStepSim.cpp 
	g++ -Ofast -I/home/anand/ -c FractionalStepSim.cpp

fractionalStepGrid.o: fractionalStepGrid.cpp 
	g++ -Ofast -I/home/anand/ -c fractionalStepGrid.cpp


fileReadingFunctions.o: fileReadingFunctions.cpp 
	g++ -Ofast -I/home/anand/ -c fileReadingFunctions.cpp

general_computation_functions.o: general_computation_functions.cpp 
	g++ -Ofast -I/home/anand/ -c general_computation_functions.cpp

grid.o: grid.cpp 
	g++ -Ofast -I/home/anand/ -c grid.cpp


testing_functions.o: testing_functions.cpp 
	g++ -Ofast -I/home/anand/ -c testing_functions.cpp

multigrid.o: multigrid.cpp 
	g++ -Ofast -I/home/anand/ -c multigrid.cpp


FracStepMultigrid.o: FracStepMultigrid.cpp 
	g++ -Ofast -I/home/anand/ -c FracStepMultigrid.cpp


clean: 
	rm -rf *o Output