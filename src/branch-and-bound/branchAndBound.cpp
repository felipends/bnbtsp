#include "branchAndBound.h"

BranchAndBound::BranchAndBound(Node* root, double** costMatrix, int dimension) {
    this->root = root;
    this->lowerBound = INFINITE;
    this->bestSolution = new Solution();
    this->costMatrix = costMatrix;
    this->dimension = dimension;
}

vector<vector<int>> BranchAndBound::getSubTours(hungarian_problem_t& p) {
    vector<vector<int>> adj(p.num_rows, vector<int>(1));
    for (int i = 0; i < p.num_rows; i++) {
        for (int j = 0; j < p.num_cols; j++) {
            if (p.assignment[i][j] == HUNGARIAN_ASSIGNED) {
                adj[i][0] = j;
            }
        }
    }

    vector<vector<int>> tours;
    vector<bool> inTuor(p.num_rows, false);
    for (int i = 0; i < p.num_rows; i++) {
        vector<int> tour;
        //se chegou no final da lista e não voltou para o i então não é um tour
        if (inTuor[i]) continue;
        tour.emplace_back(i);
        int j = i;
        while (true) {
            inTuor[j] = true;
            j = adj[j][0];
            tour.emplace_back(j);

            if (j == i) break;
        }
        tours.emplace_back(tour);
    }

    return tours;
}

vector<pair<int, int>> BranchAndBound::getForbiddenArcs(const vector<vector<int>> *subTours, int dimension, const BranchRule rule) {
    vector<pair<int, int>> forbiddenArcs;

    if (rule == SUBTOUR) {
        //find smallest subTour
        int smallestTour = 0;
        for (int i = 1; i < subTours->size(); i++) {
            if (subTours->at(i).size() < subTours->at(smallestTour).size()) {
                smallestTour = i;
            }
        }

        vector<int> chosenSubTour = subTours->at(smallestTour);
        forbiddenArcs = vector<pair<int, int>>(chosenSubTour.size() - 1);
        for (int i = 0; i < chosenSubTour.size() - 1; i++) {
            forbiddenArcs[i] = {chosenSubTour[i], chosenSubTour[i + 1]};
        }
    } else {
        //find vertex with greater degree
        vector<int> degrees(dimension, 0);
        for (auto &tour: *subTours) {
            for (int i = 0; i < tour.size(); i++) {
                degrees[tour[i]]++;
            }
        }

        int indexMaxDegree = -1;
        int maxDegree = 0;
        for (int i = 0; i < dimension; i++) {
            if (degrees[i] > maxDegree) {
                maxDegree = degrees[i];
                indexMaxDegree = i;
            }
        }

        //find arcs associated with vertex with greater degree
        forbiddenArcs = vector<pair<int, int>>();
        for (auto &tour: *subTours) {
            for (int i = 1; i < tour.size() - 1; i++) {
                if (tour[i] == indexMaxDegree) {
                    forbiddenArcs.push_back({tour[i], tour[i + 1]});
                    forbiddenArcs.push_back({tour[i-1], tour[i]});
                }
            }
        }
    }

    return forbiddenArcs;
}

Solution *BranchAndBound::solve(const BranchingStrategy strategy, const Solver solver) {
    initTree();

    double upperBound = INFINITE;

    while (!tree.empty()) {
        Node* currentNode = branching(strategy);
        solveNode(currentNode, solver);

        if (currentNode->lowerBound > upperBound) {
            delete currentNode;
            continue;
        }

        if (currentNode->feasible) {
            upperBound = min(currentNode->lowerBound, upperBound);
        } else {
            updateTree(currentNode, solver);
        }
        delete currentNode;
    }

    return bestSolution;
}

void BranchAndBound::updateTree(const Node* node, const Solver solver) {
    const auto  rule = solver == LAGRANGE ? DEGREE : SUBTOUR;
    for (auto &arc: getForbiddenArcs(&node->subTours, dimension, rule))
    {
        auto n = new Node;
        n->forbiddenArcs = node->forbiddenArcs;
        n->forbiddenArcs.emplace_back(arc);
        tree.emplace_back(n);
    }
}

void BranchAndBound::initTree() {
    for (auto &arc: root->forbiddenArcs)
    {
        auto n = new Node();
        n->forbiddenArcs = {arc};

        tree.emplace_back(n);
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

void BranchAndBound::solveHungarian(Node* node) {
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

    for (int i = 0; i < dimension; i++)
    {
        delete c[i];
    }
    delete c;

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
}

vector<vector<int>> convertTreeToTour(vector<pair<int, int>> tree) {
    vector<vector<int>> tours;
    vector<int> tour;
    for (auto &arc: tree) {
        tour.push_back(arc.first);
        tour.push_back(arc.second);
    }
    tours.push_back(tour);
    return tours;
}

void BranchAndBound::solveLagrange(Node* node) {
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

    //TODO: get lambda from node
    auto lagrangian = new Lagrangian(c, dimension, vector<double>(dimension, 0), 0);
    auto solution = lagrangian->solve();

    node->lowerBound = solution.cost;
    node->subTours = convertTreeToTour(solution.edges);
}

void BranchAndBound::solveNode(Node* node, const Solver solver) {
    switch (solver)
    {
        case LAGRANGE:
            solveLagrange(node);
            break;
        default:
            solveHungarian(node);
    }
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