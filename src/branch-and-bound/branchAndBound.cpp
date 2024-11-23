#include "branchAndBound.h"

BranchAndBound::BranchAndBound(Node* root, double** costMatrix, int dimension, double UB) {
    this->root = root;
    this->lowerBound = INFINITE;
    this->bestSolution = new Solution();
    this->costMatrix = costMatrix;
    this->dimension = dimension;
    this->UB = UB;
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

vector<pair<int, int>> BranchAndBound::getForbiddenArcs(const vector<pair<int, int>>* tree, int dimension) {
    vector<pair<int, int>> forbiddenArcs;
    //find vertex with greater degree
    vector<int> degrees(dimension, 0);
    for (auto &arc: *tree) {
        degrees[arc.first]++;
        degrees[arc.second]++;
        // cout << arc.first << " " << arc.second << endl;
    }

    int maxDegree = 0;
    int maxDegreeVertex = -1;
    for (int i = 0; i < dimension; i++) {
        if (degrees[i] > maxDegree) {
            maxDegree = degrees[i];
            maxDegreeVertex = i;
        }
    }

    // cout << maxDegree << ", " << maxDegreeVertex  << endl;
    if (maxDegree < 3) {
        return forbiddenArcs;
    }

    //extract forbidden arcs from tree
    for (auto &arc: *tree) {
        if (arc.first == maxDegreeVertex || arc.second == maxDegreeVertex) {
            forbiddenArcs.emplace_back(arc);
        }
    }

    return forbiddenArcs;
}

vector<pair<int, int>> BranchAndBound::getForbiddenArcs(const vector<vector<int>> *subTours, int dimension) {
    vector<pair<int, int>> forbiddenArcs;

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

    return forbiddenArcs;
}

void printNode(Node* node) {
    cout << "Node" << endl;
    cout << "Lower bound: " << node->lowerBound << endl;
    cout << "Feasible: " << node->feasible << endl;
    cout << "Forbidden arcs: " << endl;
    for (auto &arc: node->forbiddenArcs) {
        cout << arc.first << " " << arc.second << endl;
    }
    cout << "Sub tours: " << endl;
    for (auto &tour: node->subTours) {
        for (auto &node: tour) {
            cout << node << " ";
        }
        cout << endl;
    }

    if (node->lambda.size() > 0) {
        cout << "Lambda: " << endl;
        for (auto &l: node->lambda) {
            cout << l << " ";
        }
        cout << endl;
    }
}

Solution *BranchAndBound::solve(const BranchingStrategy strategy, const Solver solver) {
    initTree(solver);

    double upperBound = solver == LAGRANGE ? root->lowerBound : INFINITE;

    int i = 0;
    while (!tree.empty()) {
        // if (i++ == 100)
        // {
        //     break;
        // }
        Node* currentNode = branching(strategy);
        solveNode(currentNode, solver);

        //printNode(currentNode);

        if (currentNode->lowerBound > upperBound + 1) {
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
    if (solver == LAGRANGE) {
        for (auto &arc: getForbiddenArcs(&node->tree, dimension))
        {
            auto n = new Node;
            n->forbiddenArcs = node->forbiddenArcs;
            n->forbiddenArcs.emplace_back(arc);
            n->lambda = node->lambda;
            tree.emplace_back(n);
        }
    } else {
        for (auto &arc: getForbiddenArcs(&node->subTours, dimension))
        {
            auto n = new Node;
            n->forbiddenArcs = node->forbiddenArcs;
            n->forbiddenArcs.emplace_back(arc);

            tree.emplace_back(n);
        }
    }
}

void BranchAndBound::initTree(Solver solver) {
    // cout << "Root forbidden arcs: " << root->forbiddenArcs.size() << endl;
    tree.emplace_back(root);
    for (auto &arc: root->forbiddenArcs)
    {
        auto n = new Node();
        n->forbiddenArcs = vector<pair<int, int>>();//root->forbiddenArcs;
        n->forbiddenArcs.emplace_back(arc);
        n->feasible = root->feasible;
        if (solver == LAGRANGE) {
            n->lambda = root->lambda;
        }

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
        c[arc.second][arc.first] = INFINITE;
    }

    // cout << "NODE NO SOLVE LAGRANGE. ARVORE: " << tree.size() << endl;
    // printNode(node);
    // cout << "FIM DO NODE NO SOLVE LAGRANGE: " << endl;

    auto lagrangian = new Lagrangian(c, dimension, node->lambda, UB);
    auto solution = lagrangian->solve();

    // //print solution
    // cout << "Solution cost: " << solution.cost << endl;
    // for (auto &edge: solution.edges) {
    //     cout << edge.first << " " << edge.second << endl;
    // }
    // cout << "Solution printed ======================" << endl;

    node->lowerBound = solution.cost;
    node->tree = solution.edges;
    node->lambda = solution.lambda;
    node->feasible = solution.feasible;

    if (node->feasible && solution.cost < lowerBound) {
        lowerBound = solution.cost;
        bestSolution->cost = solution.cost;
        bestSolution->edges = node->tree;
    }

    for (int i = 0; i < dimension; i++)
    {
        delete c[i];
    }
    delete c;
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