# hierarchy

#### Simulation

To reproduce the simulation results, run code in `get_simu_results.R`.  Packages `ape`, `igraph`, `phangorn`, `NetworkToolbox`, and `MASS` are required.

For caterpillar trees, run
```
n_rep = 100
n = 16
h = 500
seq_length <- seq(100, 700, 100)
set.seed(1)
cater <- Get_result_binary_tree_gaussian(n_rep, n, h, seq_length, type = 'caterpillar')
```

For balanced binary trees, run
```
n_rep = 100
n = 16
h = 500
seq_length <- seq(100, 700, 100)
set.seed(1)
balance <- Get_result_binary_tree_gaussian(n_rep, n, h, seq_length, type = 'balanced')
```

For random non-binary trees, run
```
n_rep = 100
n = 16
h = 500
seq_length <- round(exp(seq(4, 10)))
set.seed(1)
random_non_binary <- Get_result_non_binary_tree_gaussian(n_rep, n, h, seq_length)
```

One could directly get simulation results by `all <- readRDS("data/all_simu_results.RDS")` and then visualize the results. To visualize the result, the `ggplot2` and `ggh4x` packages are required. Pakage `tidyr` is used to arrange the dataframe to get error bars.


#### Example of k unknown

To reproduce the example of k being unknown, run code in `robust_k_example.R`



