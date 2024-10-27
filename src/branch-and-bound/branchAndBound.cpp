#include "branchAndbound.h"

BranchAndBound::BranchAndBound(Node* root, double** costMatrix, int dimension) {
    this->root = root;
    this->lowerBound = INFINITE;
    this->bestSolution = new Solution();
    this->costMatrix = costMatrix;
    this->dimension = dimension;
}

vector<vector<int>> BranchAndBound::getSubTours(hungarian_problem_t& p) {
    vector<vector<int>> adj(p.num_rows, vector<int>());
    for (int i = 0; i < p.num_rows; i++) {
        for (int j = 0; j < p.num_cols; j++) {
            if (p.assignment[i][j] == HUNGARIAN_ASSIGNED) {
                adj[i].push_back(j);
            }
        }
    }

    vector<vector<int>> tours;
    vector<bool> inTuor(p.num_rows, false);
    for (int i = 0; i < p.num_rows; i++) {
        vector<int> tour;
        //se chegou no final da lista e não voltou para o i então não é um tour
        if (inTuor[i]) continue;
        tour.push_back(i);
        int j = i;
        while (true) {
            inTuor[j] = true;
            j = adj[j][0];
            tour.push_back(j);

            if (j == i) break;
        }
        tours.push_back(tour);
    }

    return tours;
}

vector<pair<int, int>> BranchAndBound::getForbiddenArcs(const vector<vector<int>> *subTours) {
    //find smallest subTour
    int smallestTour = 0;
    for (int i = 1; i < subTours->size(); i++) {
        if (subTours->at(i).size() < subTours->at(smallestTour).size()) {
            smallestTour = i;
        }
    }

    vector<int> chosenSubTour = subTours->at(smallestTour);
    vector<pair<int, int>> forbiddenArcs(chosenSubTour.size() - 1);
    for (int i = 0; i < chosenSubTour.size() - 1; i++) {
        forbiddenArcs[i] = {chosenSubTour[i], chosenSubTour[i + 1]};
    }

    return forbiddenArcs;
}

Solution *BranchAndBound::solve(const BranchingStrategy strategy) {
    initTree();

    double upperBound = INFINITE;

    while (!tree.empty()) {
        Node* currentNode = branching(strategy);
        currentNode = solveHungarian(currentNode);

        if (currentNode->lowerBound > upperBound) {
            delete currentNode;
            continue;
        }

        if (currentNode->feasible) {
            upperBound = min(currentNode->lowerBound, upperBound);
        } else {
            for (auto &arc: getForbiddenArcs(&currentNode->subTours))
            {
                auto n = new Node;
                n->forbiddenArcs = currentNode->forbiddenArcs;
                n->forbiddenArcs.push_back(arc);
                tree.push_back(n);
            }
        }
        delete currentNode;
    }

    return bestSolution;
}

void BranchAndBound::initTree() {
    for (auto &arc: root->forbiddenArcs)
    {
        auto n = new Node();
        n->forbiddenArcs = {arc};

        tree.push_back(n);
    }
}

void print_tours(vector<vector<int>> &tours) {
    cout << "Tours" << endl;
    for (auto &tour : tours) {
        for (auto &node : tour) {
            cout << node << " ";
        }
        cout << endl;
    }
}

Node* BranchAndBound::solveHungarian(Node* node) {
    double** c = new double*[dimension];
    for (int i = 0; i < dimension; i++) {
        c[i] = new double[dimension];
        for (int j = 0; j < dimension; j++) {
            c[i][j] = costMatrix[i][j];
        }
    }

    for (auto &arc: node->forbiddenArcs) {
        c[arc.first][arc.second] = INFINITE;
    }

    hungarian_problem_t p;
    int mode = HUNGARIAN_MODE_MINIMIZE_COST;
    hungarian_init(&p, c, dimension, dimension, mode);

    double obj_value = hungarian_solve(&p);

    node->lowerBound = obj_value;
    node->subTours = getSubTours(p);

    if (node->subTours.size() == 1) {
        node->feasible = true;
        if (obj_value < lowerBound) {
            print_tours(node->subTours);
            lowerBound = obj_value;
            bestSolution->cost = obj_value;
            bestSolution->tour = node->subTours[0];
        }
    } else {
        node->feasible = false;
    }

    hungarian_free(&p);
    delete [] c;
    return node;
}

Node *BranchAndBound::branching(const BranchingStrategy strategy) {
    Node* currentNode = nullptr;
    switch (strategy) {
    case BFS:
        currentNode = tree.front();
        tree.pop_front();
        break;
    default:
        currentNode = tree.back();
        tree.pop_back();
        break;
    }

    return currentNode;
}