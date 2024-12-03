#include <iostream>

#include "branch-and-bound/branchAndBound.h"
using namespace std;

#include "instance-reader/Data.h"
#include "hungarian/Hungarian.h"

int main(int argc, char** argv) {

	Data * data = new Data(argc, argv[1]);
	data->readData();

	double **cost = new double*[data->getDimension()];
	for (int i = 0; i < data->getDimension(); i++){
		cost[i] = new double[data->getDimension()];
		for (int j = 0; j < data->getDimension(); j++){
			cost[i][j] = data->getDistance(i,j);
		}
	}

	vector<Solver> solvers = {HUNGARIAN, LAGRANGE};
	vector<BranchingStrategy> strategies = {DFS, BFS, BEST_BOUND};

	auto bnb = new BranchAndBound(cost, data->getDimension(), atof(argv[4]));
	// cout << "CRIOU" << endl;
	auto solution = bnb->solve(strategies[atoi(argv[2])], solvers[atoi(argv[3])]);
	// cout << "RODOU O SOLVE" << endl;

	cout << "Solution: " << solution->cost << endl;

	for (auto s: solution->tour) {
		cout << s << " ";
	}
	cout << endl;

	for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	delete [] cost;
	delete data;

	return 0;
}