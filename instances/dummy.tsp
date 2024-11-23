NAME: dummy_instance_extended
TYPE: TSP
COMMENT: Exemplo de instância TSP com matriz triangular inferior estendida com um nó adicional (todos valores 0).
DIMENSION: 5
EDGE_WEIGHT_TYPE: EXPLICIT
EDGE_WEIGHT_FORMAT: UPPER_DIAG_ROW
EDGE_WEIGHT_SECTION
0 30 26 50 40
   0 24 40 50
      0 24 26
         0 30
            0
EOF