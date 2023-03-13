all: Output

Output: main.o FractionalStepSim.o fileReadingFunctions.o fractionalStepGrid.o  general_computation_functions.o grid.o testing_functions.o multigrid.o FracStepMultigrid.o 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 main.o FractionalStepSim.o fileReadingFunctions.o fractionalStepGrid.o  general_computation_functions.o grid.o testing_functions.o multigrid.o FracStepMultigrid.o -o Output

main.o: main.cpp
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c main.cpp 

FractionalStepSim.o: FractionalStepSim.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c FractionalStepSim.cpp

fractionalStepGrid.o: fractionalStepGrid.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c fractionalStepGrid.cpp


fileReadingFunctions.o: fileReadingFunctions.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c fileReadingFunctions.cpp

general_computation_functions.o: general_computation_functions.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c general_computation_functions.cpp

grid.o: grid.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c grid.cpp


testing_functions.o: testing_functions.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c testing_functions.cpp

multigrid.o: multigrid.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c multigrid.cpp


FracStepMultigrid.o: FracStepMultigrid.cpp 
	g++ -Ofast -std=c++17 -I/home/anand/ -I/opt/homebrew/include/eigen3 -c FracStepMultigrid.cpp


clean: 
	rm -rf *o Output