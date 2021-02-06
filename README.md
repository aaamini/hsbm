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
|method_name | aggregate_nmi| slicewise_nmi| elapsed_time|
|:-----------|-------------:|-------------:|------------:|
|HSBM        |        0.9766|        0.9804|        0.477|
|DP-SBM      |        0.1831|        0.8920|        1.450|
|SC-sliced   |        0.0315|        0.2006|        0.089|
|SC-avg      |        0.0018|        0.0049|        0.034|
```
Here, `DP-SBM` refers to the algorithm that runs the Dirichlet Process SBM separately on each layer. `SC-sliced` runs the spectral clustering separately on each layer (i.e., slice).  `SC-avg` run the spectral clustering on the average of the adjacency matrices from all layers. 

Changelog:
- 2/5/2021: The code was rewritten from scratch resulting in a much faster and more stable sampler. The code is optimized to work with sparse networks and runs ~ 100x fast than the old code (now moved to `old_code` folder). The sampler also mixes fast, on average in about 10-50 iterations. The new code is now provided in the R package `hsbm` available under the sub-folder `hsbm-package`. 


