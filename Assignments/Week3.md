**Goal:** The goal is to become familiar with the compute cluster at WCM.

1. scp Week3.R to aphrodite  &&
2. create / scp Week3.sh to aphrodite  &&
3. Install all needed packages in the server (using qlogin) -- only the first time!  &&
4. qsub Week3.sh   &&
5. scp results back to your laptop
6. open results with Rstudio to make sure it worked.


**Goal:** The goal is to update the C++ code for mini-batch clustering in the clusterR package to pull random subsets from hdf5 files instead of a matrix in memory.

Basic outline of steps:
1) Learn how to use C++ commands to read a hdf5 file. Like the very first assignment, create a small function that takes as input a hdf5 file name, divides the data into chunks, and returns the means of specific chunks

2) Familiarize your self with the idea of mini-Batch clustering (mainly via a quick tutorial from meeting with Davide)

3) Identify the C++ code in the clusterR package (get the "source" code of the package if there's not a github repository) and make a copy of it to play with in the beachball package. Make sure the code is in the right place for the package structure. There may be more than one C++ function necessary. Make sure that there is a R function in beachball that can successfully call the miniBatch C++ code (i.e. independently of clusterR package) and that it gives the same answer.

4) Change the miniBatch C++ code to modularize the subsetting of the data. e.g. make a function f that takes indices and full data matrix as input and returns the necessary data matrix (exact format should be what makes most sense with existing C++ code). 


5) Now, create a alternative C++ module that takes the same input as the function f in your mini-batch code, only now takes the input of full data as instead a hdf5 file name. This should all be using C++ commands. 

6) Integrate this hdf5 version of the f function into the mini-batch code. It may be that the best way to do this is to write two input functions, one for hdf5MiniBatch and one for matrixMiniBatch, and then the two functions call different functions f for getting the subset data, but then both call the *same* function to do the clustering part of the function. But DON'T COPY CODE. If you need the same code twice, it needs to be made a function that gets called (and of course with the fewest number of functions possible).
7) (if needed) Update the R code that calls the C++ code so that it appropriately switches between the hdf5 version and the in-memory version of the C++ code 