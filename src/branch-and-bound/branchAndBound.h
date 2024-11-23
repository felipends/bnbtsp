#ifndef BRANCH_AND_BOUND_H
#define BRANCH_AND_BOUND_H
#include <list>
#include <iostream>

#include "../hungarian/Hungarian.h"
#include "../lagrangian-relaxation/lagrangian.h"
#define INFINITE 99999999

using namespace std;

void print_tours(vector<vector<int>> &tours);

typedef struct
{
    vector<pair<int, int>> forbiddenArcs;
    vector<pair<int, int>> tree;
    vector<vector<int>> subTours;
    double lowerBound;
    int chosenSubTour;
    bool feasible;

    //lagrange
    vector<double> lambda;
} Node;

typedef struct
{
    vector<int> tour;
    vector<pair<int, int>> edges;
    double cost;
} Solution;

enum BranchingStrategy
{
    DFS,
    BFS,
    BEST_BOUND
};

enum Solver
{
    HUNGARIAN,
    LAGRANGE
};

enum BranchRule
{
    DEGREE,
    SUBTOUR
};

class BranchAndBound {
public:
    BranchAndBound(Node* root, double** costMatrix, int dimension, double UB = INFINITE);

    Solution* solve(BranchingStrategy strategy = DFS, Solver solver = HUNGARIAN);

    static vector<vector<int>> getSubTours(hungarian_problem_t &p);
    static vector<pair<int, int>> getForbiddenArcs(const vector<vector<int>> *subTours, int dimension);
    static vector<pair<int, int>> getForbiddenArcs(const vector<pair<int, int>> *tree, int dimension);
private:
    list<Node*> tree;
    Node* root;
    double lowerBound;
    Solution* bestSolution;
    int dimension;
    double** costMatrix;
    double UB;

    void initTree(Solver solver);
    Node* branching(BranchingStrategy strategy = DFS);

    void solveHungarian(Node* node);
    void solveLagrange(Node* node);
    void solveNode(Node* node, Solver solver);
    void updateTree(const Node* node, Solver solver);
};

#endif