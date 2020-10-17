# Parallel-Computing-Assignment
Parallel Computing Assignment Shortest Path Algorithm Using Matrix Multiplication

Assumption : Number of nodes of graph is divisible by number of process. 

An undirected graph is created randomly for each run, with edge weight 1 and is written into a file graph.txt. If graph size is small(<=32) the graph is displayed on the terminal otherwise not.

The shortest path algorithm is run on the graph first parallely with n process than on a single process to verify the output and compare the time of execution.

Observation : For small graph serial algorithm gives better performance, the reason is commmunication overhead in the parallel algorithm, but for large graphs parallel algorithm gives far better performance then serial algorithm taking much less time to execute the program and produce the output compared to serial algorithm.

