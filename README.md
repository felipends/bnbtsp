# Branch and Bound para resolver o TSP

Implementação inclui a opção de resolver o TSP com a estratégia de branching BFS, DFS e BEST BOUND.

## Execução
```bash
cmake -B build
make -C build
./build/tsp <arquivo de entrada> <branching strategy (0 - DFS|1 - BFS)> <solver (0 - HUNGARIAN|1 - LAGRANGIAN)> <UB - number: upperbound conhecido para instância>
```