# MATLAB Toolbox for fdipa, feasible direction interior-point algorithm

MATLAB library to find minimum of nonlinear second-order cone programming problem, using the feasible direction interior-point algorithm (fdipa) as presented in [1].

## Download 
Option 1: From a MATLAB command prompt run:
```
git clone https://github.com/fdipaSOC/NSOCP
```
Option 2: run script "fdipa_setup.m" included in the library folder

## Install in MATLAB
From a MATLAB command window:

```
addpath( genpath('AbsolutePathToToolbox') )
```
where `AbsolutePathToToolbox` is the name of the absolute path where you cloned this toolbox. 
Adding this command to your MATLAB `startup.m` file will make sure these tools are avalible every time you
use MATLAB.

## Documentation

[Documentation PDF](documentation/doc_fdipa.pdf)

## Contents
The scripts are organized in the following folders:

* `fdipa`: source files for fdipa algorithm. 
  * `fdipa/fdipa.m`: solves a random linear programming problem using either interface.
  * `fdipa/fdipaQuad.m`: solves a random linear programming problem using either interface.
* `documentation`: source files for the documentation
* `example1_quadratic`: example of application of fdipaQuad to a quadratic problem with linear cone constraints.
* `example2_kato`: Kato-Fukushima example of for nonlinear second-order cone programs as presented in [2].
* `example3_miao`: Miao-Chen-Ko example of nonlinear convex programs with second-order cone constraints as presented in [3].
* `example4_hyf_nonlinear`: Example of nonlinear problem with linear constraints used in Example 4.2 of [4] (also  used in Ex. 5.3 of [3]).
* `example5_dist_ellip`: application of fdipaQuad to the problem of distance between two ellipses.
* `example6_cubic`: example of cubic objective function with linear constraints.
* `example7_rotated_conic`: example of a quadratic objective with rotated conic constraints.
* `example8_svm`: Application to support vector machines (Requires CVX).

## Citing
This repository has been developed as part of the following article (citation given below). We would appreciate it if you would please cite the following paper if you found the library useful for your work.

```
@article{fdipa2024,
    title = {FDIPA-SOC: A Matlab package for nonlinear Second-Order Cone programs},
    author = {Canelas, A. and Carrasco, M and L\'opez, J and Paduro, E.},
    journal = {2024},
    volume = {},
    pages = {},
    year = {},
    publisher={},
    doi = {},
    url = {}
}
```

## References


[1] Alfredo Canelas, Miguel Carrasco & Julio Lopez (2019) A feasible direction algorithm for nonlinear second-order cone programs, *Optimization Methods and Software*, 34:6, 1322-1341, 
[https://doi.org/10.1080/10556788.2018.1506452](https://doi.org/10.1080/10556788.2018.1506452)

[2] Kato, H., Fukushima, M. An SQP-type algorithm for nonlinear second-order cone programs. *Optimization Letters 1*, 129–144 (2007). 
[https://doi.org/10.1007/s11590-006-0009-2](https://doi.org/10.1007/s11590-006-0009-2)

[3] Xinhe Miao, Jein-Shan Chen, Chun-Hsu Ko. A smoothed NR neural network for solving nonlinear convex programs with second-order cone constraints, *Information Sciences*, Volume 268, 2014, p 255-270, 
[https://doi.org/10.1016/j.ins.2013.10.017](https://doi.org/10.1016/j.ins.2013.10.017)

[4] C. Kanzow, I. Ferenczi, and M. Fukushima. On the local convergence of semismooth newton methods for linear and nonlinear second-order cone programs without  strict complementarity. *SIAM J. Optim.*, 20(1):297-320, 2009.
[https://doi.org/10.1137/060657662](https://doi.org/10.1137/060657662)

## License

[GPLv3](LICENSE)