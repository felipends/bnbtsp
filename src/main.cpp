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
	Node* root = new Node();
	if (solvers[atoi(argv[3])] == HUNGARIAN)
	{
		hungarian_problem_t p;
		int mode = HUNGARIAN_MODE_MINIMIZE_COST;
		hungarian_init(&p, cost, data->getDimension(), data->getDimension(), mode); // Carregando o problema

		double obj_value = hungarian_solve(&p);
		cout << "Obj. value: " << obj_value << endl;

		vector<vector<int>> tours = BranchAndBound::getSubTours(p);
		print_tours(tours);

		root->lowerBound = obj_value;
		root->feasible = false;
		root->subTours = tours;
		root->forbiddenArcs = BranchAndBound::getForbiddenArcs(&tours, data->getDimension());
		hungarian_free(&p);
	} else if (solvers[atoi(argv[3])] == LAGRANGE) {
		auto lambda = vector<double>(data->getDimension(), 0);
		auto forbiddenArcs = vector<pair<int, int>>();
		auto lagrangian = new Lagrangian(cost, forbiddenArcs, data->getDimension(), lambda, atof(argv[4]));
		auto solution = lagrangian->solve();
		cout << "LB Root: " << solution.cost << endl;
		//print tree
		for (auto &edge: solution.edges) {
			cout << edge.first << " " << edge.second << endl;
		}

		root->lowerBound = solution.cost;
		root->feasible = solution.feasible;
		root->tree = solution.edges;
		root->forbiddenArcs = BranchAndBound::getForbiddenArcs(&root->tree, data->getDimension());
		root->lambda = solution.lambda;
	}

	auto bnb = new BranchAndBound(root, cost, data->getDimension(), atof(argv[4]));
	vector<BranchingStrategy> strategies = {DFS, BFS, BEST_BOUND};
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