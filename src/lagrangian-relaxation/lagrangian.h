#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H
#include <vector>

#define E_MIN 0.000001
#define K_MAX 30

#define INFINITE 999999999

using namespace std;

typedef struct
{
    vector<pair<int, int>> edges;
    vector<double> lambda;
    double cost;
    bool feasible;
} LagrangianSolution;

class Lagrangian {
public:
    Lagrangian(double** costMatrix, vector<pair<int, int>>& forbiddenArcs, int dimension, vector<double> lambda, double UB);

    LagrangianSolution solve() const;
    bool isFeasible(const vector<int>& verticesDegrees) const;
    double getSubgradient(const vector<int>& verticesDegrees) const;
    static vector<int> getVerticesDegrees(const vector<pair<int, int>>& edges, int dimension);

private:
    vector<vector<double>> costMatrix;
    static vector<vector<double>> updateCostMatrix(vector<vector<double>> costMatrix, const vector<double>& lambda);
    int dimension;
    vector<double> lambda;
    double UB;
    vector<pair<int, int>> forbiddenArcs;
    void penalizeForbiddenArcs(vector<vector<double>>& costMatrix) const;
};

#endif //LAGRANGIAN_H
