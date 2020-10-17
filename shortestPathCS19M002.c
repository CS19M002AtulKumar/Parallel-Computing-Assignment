/* Shortest path using matrix multiplication, submitted by CS19M002 */
// ASSUMPTION : NUMBER OF NODES IS DIVISIBLE BY NUMBER OF PROCESS.
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 512		  			// Number of nodes in the graph
#define SOURCE 0					// Source node
#define DESTINATION SIZE-1			// Destination node


/* parallelGraph to run the parallel algorithm, serialGraph to run the shortest path algorithm serially to verify the output. */
int adjParallel[SIZE][SIZE], adjSerial[SIZE][SIZE], result[SIZE][SIZE];


/* Writes the graph into a file */
void writeGraphIntoFile(int graph[SIZE][SIZE]){
	FILE *filePointer;
	filePointer = fopen("graph.txt", "w");
	for(int i=0; i<SIZE; i++){
		for(int j=0; j<SIZE; j++){
			fprintf(filePointer,"%d ", graph[i][j]);
		}
		fprintf(filePointer, "\n");
	}
	fclose(filePointer);
}


/* Creates a random undirected graph with edge weight one. */
void createGraph(int parallelGraph[SIZE][SIZE], int serialGraph[SIZE][SIZE])
{
  srand(time(0));
  int val;
  for(int i=0; i<SIZE; i++){
  	for(int j=0; j<i; j++){
  		if((rand()%100) < (rand()%30)){
  			val = 1;
  		}
  		else{
  			val = 0;
  		}  		
  		parallelGraph[i][j] = val;
  		parallelGraph[j][i] = val;
  		serialGraph[i][j] = val;
  		serialGraph[j][i] = val;
  	}
  }
}


/* print the graph on the terminal */
void printGraph(int arr[SIZE][SIZE])
{
  for (int i=0; i<SIZE; i++) {
    for (int j=0; j<SIZE; j++){
      printf("%d ", arr[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void copyAdj(int arr[SIZE][SIZE], int res[SIZE][SIZE]){
	for(int i=0; i<SIZE; i++){
    	for(int j=0; j<SIZE; j++){
    		arr[i][j] = res[i][j];
    	}
    }
}

int main(int argc, char *argv[])
{
  int rank, noOfProcess, from, to, length=1;
  int flag = 0;
  
  MPI_Init (&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
  MPI_Comm_size(MPI_COMM_WORLD, &noOfProcess); 


  if (SIZE%noOfProcess!=0) {
    if (rank==0){
    	printf("Number of nodes should be divisible by number of process.\n");
    } 
    MPI_Finalize();
    exit(0);
  }
  
  if(SOURCE>=SIZE || DESTINATION>=SIZE){
  	if(rank==0)
  	printf("Source node and Destination node should be between (0,%d)\n", SIZE-1);
  	MPI_Finalize();
  	exit(0);
  }
  
  if (rank==0) {
    createGraph(adjParallel, adjSerial);
    // Write the graph into the file
  	writeGraphIntoFile(adjParallel);
  	// Print the graph on terminal if the graph size is small
  	if(SIZE <= 32 ){
  		printGraph(adjParallel);
  	}
  	else{
  		printf("Graph written into file graph.txt. Too large to print on terminal.\n");
  	}
  	printf("Source node : %d\n", SOURCE);
  	printf("Destination node : %d\n", DESTINATION);
  	printf("Number of process : %d\n", noOfProcess);
  	// implementing shortest path algorithm with parallel computation.
    if(adjParallel[SOURCE][DESTINATION]!=0){
    	flag = 1;
    }
  }
  
  clock_t beginParallel = clock();  
  if(flag==0){
	  from = rank * SIZE/noOfProcess;
	  to = (rank+1) * SIZE/noOfProcess;
  }
  
  
  MPI_Bcast (&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  while(flag==0 && length!=SIZE){
  
  	  length++;
	  MPI_Bcast (adjParallel, SIZE*SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	  
	  for (int i=from; i<to; i++) 
		for (int j=0; j<SIZE; j++) {
		  result[i][j]=0;
		  for (int k=0; k<SIZE; k++)
		result[i][j] += adjParallel[i][k]*adjParallel[k][j];
		}

	  MPI_Gather (result[from], SIZE*SIZE/noOfProcess, MPI_INT, result, SIZE*SIZE/noOfProcess, MPI_INT, 0, MPI_COMM_WORLD);

	  if (rank==0) {
		copyAdj(adjParallel, result);
		if(adjParallel[SOURCE][DESTINATION]!=0){
			flag = 1;
		}
	  }
	  
	  MPI_Bcast (&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);	  
	  
  }
  
  if(rank==0){
  		clock_t endParallel = clock();
		double timeSpentParallel = (double)(endParallel - beginParallel) / CLOCKS_PER_SEC;
  		if(length==SIZE){
  			printf("Implementing parallel algorithm path : no Path From Source To Destination!!!!\n");
  		}
  		else{
  			printf("Implementing parallel algorithm path length : %d\n", length);
  		}
  		printf("Time spent in parallel algorithm execution with %d processes : %lf\n",noOfProcess, timeSpentParallel );
  		
  		// Now implementing shortest path algorithm in a single process to compare time and verify result.
  		
  		int pathLength = 1;
  		clock_t beginSerial = clock();
  		while(adjSerial[SOURCE][DESTINATION]==0 && pathLength!=SIZE){
  			for(int i=0; i<SIZE; i++){
  				for(int j=0; j<SIZE; j++){
  					result[i][j]=0;
  					for(int k=0; k<SIZE; k++){
  						result[i][j] += adjSerial[i][k]*adjSerial[k][j];
  					}
  				}
  			}
  			copyAdj(adjSerial, result);
  			pathLength++;
  			
  		}
  		clock_t endSerial = clock();
		double timeSpentSerial = (double)(endSerial - beginSerial) / CLOCKS_PER_SEC;
  		if(pathLength==SIZE){
  			printf("Implementing serial algorithm path : no Path From Source To Destination!!!!\n");
  		}
  		else{
  			printf("Implementing serial algorithm path length : %d\n", pathLength);
  		}
  		printf("Time spent in serial algorithm execution : %lf\n", timeSpentSerial );
  }

  MPI_Finalize();

  return 0;
}




