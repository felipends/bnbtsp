#ifndef BRANCH_AND_BOUND_H
#define BRANCH_AND_BOUND_H
#include <list>
#include <iostream>

#include "../hungarian/Hungarian.h"
#define INFINITE 99999999

using namespace std;

void print_tours(vector<vector<int>> &tours);

typedef struct
{
    vector<pair<int, int>> forbiddenArcs;
    vector<vector<int>> subTours;
    double lowerBound;
    int chosenSubTour;
    bool feasible;
} Node;

typedef struct
{
    vector<int> tour;
    double cost;
} Solution;

enum BranchingStrategy
{
    DFS,
    BFS,
    BEST_BOUND
};

class BranchAndBound {
public:
    BranchAndBound(Node* root, double** costMatrix, int dimension);

    Solution* solve(BranchingStrategy strategy = DFS);

    static vector<vector<int>> getSubTours(hungarian_problem_t &p);
    static vector<pair<int, int>> getForbiddenArcs(const vector<vector<int>> *subTours);
private:
    list<Node*> tree;
    Node* root;
    double lowerBound;
    Solution* bestSolution;
    int dimension;
    double** costMatrix;

    void initTree();
    Node* branching(BranchingStrategy strategy = DFS);

    void solveHungarian(Node* node);
};

#endif