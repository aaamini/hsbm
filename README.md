# hsbm
*Hierarchical stochastic block model (HSBM)*

The code implements the algorithm discussed in [Hierarchical Stochastic Block Model for Community Detection in Multiplex Networks](https://arxiv.org/abs/1904.05330).

The code now relies on the [`nett` package](https://github.com/aaamini/nett) for some of its functionality. 

How to use:
- Install the `nett` package from the above link.
- Install the `hsbm` package from this repository by issuing the following command:
```
devtools::install_github("aaamini/hsbm", subdir = "hsbm_package")
```
- Run the `benchmark.R` in the root of the repository. The ouput would be something like this:
```
|method    | aggregate_nmi| slicewise_nmi| elapsed_time|
|:---------|-------------:|-------------:|------------:|
|HSBM      |        0.9766|        0.9804|        0.238|
|DP-SBM    |        0.1831|        0.8920|        1.197|
|SC-sliced |        0.0315|        0.2006|        0.098|
|SC-avg    |        0.0018|        0.0049|        0.035|
|SC-ba     |        0.0029|        0.0765|        0.086|
|SC-omni   |        0.0038|        0.0434|        0.199|
|PisCES    |        0.1016|        0.1408|        0.125|
|PisCES-sh |        0.0347|        0.1406|        0.094|
```
Besides HSBM, the following methods are implemented:

- `DP-SBM` refers to the algorithm that runs the Dirichlet Process SBM separately on each layer. 
- `SC-sliced` runs the spectral clustering separately on each layer (i.e., slice).  
- `SC-avg` run the spectral clustering on the average of the adjacency matrices from all layers. 
- `SC-bc` is the [biased-adjusted spectral clustering](https://arxiv.org/abs/2003.08222). 
- `SC-omni` is the spectral clustering based on the [omnibus embedding](https://arxiv.org/abs/1705.09355) (with minor modification)
- `PisCES` The [PisCES](https://www.pnas.org/content/115/5/927) algorithm that solves an optimization problem that smooths out spectral projection matrices across time.
- `PisCES-sh` Our variation on the PisCES with shared k-means initialization. See the HSBM paper (updated version).

Changelog:
- 2/5/2021: The code was rewritten from scratch resulting in a much faster and more stable sampler. The code is optimized to work with sparse networks and runs ~ 100x faster than the old code (now moved to `old_code/` folder). The sampler also mixes fast, on average in about 10 to 50 iterations. The new code is now provided in the R package `hsbm` available under the sub-folder `hsbm-package`. 


