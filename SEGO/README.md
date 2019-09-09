### Outline

- [Compile](#Compile)
- [Usage](#Usage)
  - [For example](#For-example)

### Compile
Install GCC 4.7 or newer version in the [link](https://sourceforge.net/projects/mingw-w64/files/). Execute the command ```make``` in the path of the Makefile to compile the code.

    make

### Usage
Ten parameters are required for the program. Four parameters are for global optimization problem. 
1. Output path
2. Number of runs
3. Dataset number
4. Number of dimensions

Six parameters are for SEGO.
5. Searcher rate 
6. Number of regions in the initialization
7. Goods rate
8. Number of memories
9. Tournament rate 
10. Reduction period

#### For Example
Output path is **result**, number of runs is 51. Solve the test function 5 and *D*=30. It will create a file named **se_5_30.txt** in the **result/** folder. The parameter of the searcher rate is 21, number of regions in the initialization is 4, goods rate is 3.5, number of memories is 7 tournament rate is 0.011, and reduction period is 0.45.

    se.exe result 51 5 30 21 4 3.5 7 0.0011 0.45
