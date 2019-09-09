### Search economics, a novel metaheuristic algorithm, for the global optimization problem
This code implements search economics (SE) to solve the global optimization problem.

- [Introduction](#Introduction)
- [Folder](#Folder)

### Introduction
[Search economics (SE)](https://doi.org/10.1109/SMC.2015.447) is a novel metaheuristic algorithm proposed in 2015. The two main differences between SE and traditional metaheuristic are dividing the search space and investing computation resource depending on the potential of each subspace. A modified SE for solving the global optimization problem, called search economics for global optimization (SEGO).  As a pioneered approach for the global optimization problem, the modified SE includes three new mechanisms---parameter adaptation, linear population size reduction, and new search strategy---and two modified mechanisms---the way the search space is divided and the way the expected value is computed. Compare it with four state-of-the-art search algorithms (i.e., [jSO](https://doi.org/10.1109/CEC.2017.7969456), [EBOwithCMAR](https://doi.org/10.1109/CEC.2017.7969524), [ELSHADE-SPACMA](https://doi.org/10.13140/RG.2.2.33283.20005), and [L-SHADE-RSP](https://doi.org/10.1109/CEC.2018.8477977)) for solving 29 functions from CEC2017 in the [thesis](http://www.airitilibrary.com/Publication/alDetailedMesh1?DocID=U0005-0408201921375500).

### Folder
There are three folders here. **CEC2017** is the benchmark code and document of the CEC competitions in 2017. The benchmark is able to be downloaded in the [link](https://www.ntu.edu.sg/home/EPNSugan/index_files/CEC2017/CEC2017.htm). **SEGO** is the code of SEGO.