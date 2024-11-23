#include "lagrangian.h"

#include "mst/Kruskal.h"

Lagrangian::Lagrangian(double** costMatrix, int dimension, vector<double> lambda, double UB) {
    this->UB = UB;
    this->lambda = lambda;
    this->costMatrix = vector<vector<double>>();
    for (int i = 0; i < dimension; i++) {
        vector<double> costVector(dimension);
        for (int j = 0; j < dimension; j++) {
            costVector[j] = costMatrix[i][j];
        }
        this->costMatrix.push_back(costVector);
    }

    this->dimension = dimension;
}

void printEdges(const vector<pair<int, int>>& edges, vector<double> lambda) {
    cout << "Edges:" << endl;
    for (auto& edge: edges) {
        cout << edge.first + 1 << " " << edge.second + 1 << endl;
    }

    cout << "Lambda: " << endl;
    for (auto& l: lambda) {
        cout << l << " ";
    }
    cout << endl;
}

void printMatrix(vector<vector<double>> matrix) {
    cout << "Matrix:" << endl;
    for (auto& row: matrix) {
        for (auto& cell: row) {
            cout << cell << " ";
        }
        cout << endl;
    }
}

LagrangianSolution Lagrangian::solve() const
{
    auto lambda = this->lambda;
    double epsilon = 1;
    int k = 0;

    LagrangianSolution bestSolution{
        .edges = vector<pair<int, int>>(),
        .lambda = vector<double>(dimension, 0),
        .cost = 0
    };

    while(true)
    {
        LagrangianSolution currentSolution{
            .edges = vector<pair<int, int>>(),
            .lambda = lambda,
            .cost = 0
        };
        const auto penalizedMatrix = updateCostMatrix(costMatrix, lambda);
        Kruskal kruskal(penalizedMatrix);
        currentSolution.cost = kruskal.MST(dimension);
        currentSolution.edges = kruskal.getEdges();

        // printMatrix(penalizedMatrix);
        // printEdges(currentSolution.edges, currentSolution.lambda);
        // cout << "Cost: " << currentSolution.cost << endl;

        if (isFeasible(getVerticesDegrees(currentSolution.edges, dimension))) {
            currentSolution.feasible = true;
        }

        if (currentSolution.cost > bestSolution.cost) {
            bestSolution = currentSolution;
            k = 0;
        } else {
            k++;
            if (k >= K_MAX) {
                k = 0;
                epsilon /= 2;
            }
        }

        auto degrees = getVerticesDegrees(currentSolution.edges, dimension);
        const double step = epsilon * (UB - currentSolution.cost)/getSubgradient(degrees);
        for (int i = 1; i < dimension; i++) {
            lambda[i] += step * (2 - degrees[i]);
        }
        if (epsilon < E_MIN || currentSolution.cost >= UB - numeric_limits<double>::epsilon() || currentSolution.feasible) break;
    }

    return bestSolution;
}

bool Lagrangian::isFeasible(const vector<int>& verticesDegrees) const {
    for (int i = 0; i < dimension; i++) {
        if (verticesDegrees[i] > 2) return false;
    }

    return true;
}

double Lagrangian::getSubgradient(const vector<int>& verticesDegrees) const {
    double subgradient = 0;
    for (int i = 0; i < dimension; i++) {
        subgradient += (2 - verticesDegrees[i]) * (2 - verticesDegrees[i]);
    }

    return subgradient;
}

vector<int> Lagrangian::getVerticesDegrees(const vector<pair<int, int>>& edges, int dimension) {
    auto degrees = vector<int>(dimension, 0);
    for (auto& edge: edges) {
        degrees[edge.first]++;
        degrees[edge.second]++;
    }

    return degrees;
}

vector<vector<double>> Lagrangian::updateCostMatrix(vector<vector<double>> costMatrix, const vector<double>& lambda) {
    for (int i = 0; i < costMatrix.size(); i++) {
        for (int j = 0; j < costMatrix[i].size(); j++) {
            if (i != j) costMatrix[i][j] -= lambda[i] + lambda[j];
            else costMatrix[i][j] = 0;
        }
    }

    return costMatrix;
}