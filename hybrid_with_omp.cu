#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cudpp.h>
#include<limits.h>
#include <sys/time.h>
#include <omp.h>

#define NO_OF_THREADS_PER_BLOCK  1024
#define OMP_NUM_THREADS	16

float f;
unsigned int  noOfEdges;
unsigned int  noOfVertices;

unsigned int *vertices;
unsigned int *edges;
unsigned int *weights;
unsigned int *d_size;
unsigned int *d_edgeListSize;
unsigned int *d_vertexListSize;
unsigned int *segmentedMinScanInput;
unsigned int *d_segmentedMinScanInput;
unsigned int *d_segmentedMinScanOutput;
unsigned int *d_previousIDs;
unsigned int *d_successorArray;
unsigned int *d_successorArrayTemp;
unsigned int *d_indices;
unsigned int *d_edgeMap;
unsigned int *d_edgeMapCopy;
unsigned int *d_edgesCopy;
unsigned int *d_edgeIndices;
unsigned int *d_superVertexID;
unsigned int *d_superEdgeId;
unsigned int *d_MSTOutput;
unsigned int *h_MSTOutput;

unsigned int *d_edges;
unsigned int *d_vertices;
unsigned int *d_weights;
unsigned int *d_edgeFlagArray;
unsigned int *d_vertexFlagArray;

unsigned int  noOfEdgesOriginal;
unsigned int  noOfVerticesOriginal;
int *d_pickArray;

CUDPPHandle theCudpp;
CUDPPHandle             segmentedScanPlan_min;
CUDPPConfiguration	segmented_min_scan_config;

CUDPPHandle 		scanPlan;
CUDPPConfiguration 	scan_config;

CUDPPHandle 		sortPlan;
CUDPPConfiguration 	config_sort;

/* Append vertexid and edge into a single integer of an array*/
__global__ void mergeEdgeAndWeight(unsigned int *d_segmentedMinScanInput, unsigned int *d_vertices, unsigned int *d_weight, unsigned int *d_edges, unsigned int noOfEdges) 
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfEdges) {
		unsigned int temp = d_weight[index];
		d_segmentedMinScanInput[index] = (temp<<22) | d_edges[index];
	}
}

/* initialise all entries of array pointed by d_array of given size to 0*/
__global__ void initArray(unsigned int *d_Array, unsigned int size) 
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	if(index < size) {
		d_Array[index] = 0;
	}
}

__global__ void initArray1(unsigned int *d_Array, unsigned int size, int t)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
        if(index < size && index >= t)
                d_Array[index] = 0;
}

__global__ void printArr(unsigned int *d_arr, unsigned int size)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	
	if (index < size) {
		printf("%d ", d_arr[index]);
	}
	printf("\n");
}

/* creates a flag array for segmented scan. Sets to 1 the index from where outgoing vertex starts*/
__global__ void markSegment(unsigned int *d_edgeFlagArray, unsigned int *d_vertex, unsigned int *d_edges, unsigned int noOfVertices) 
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	if(index < noOfVertices) {
		d_edgeFlagArray[d_vertex[index]] = 1;
	}
}

/*prints new edge and vertex size*/
__global__ void print(unsigned int *d_edgeListSize, unsigned int *d_vertexListSize)
{
	printf("Edges: %d, Vertices %d \n", *d_edgeListSize, *d_vertexListSize);

}

/*creates successor array*/
__global__ void createSuccArray(unsigned int *d_successorArray, unsigned int *d_vertices, unsigned int *d_segmentedMinScanInput, unsigned int *d_segmentedMinScanOutput, unsigned int noOfVertices, unsigned int noOfEdges)
{
        unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	unsigned int minEdgeIndex;

	if(index < noOfVertices) {
		//index is same as vertex ID
		if (index == noOfVertices-1)
			minEdgeIndex = noOfEdges - 1;
		else
			minEdgeIndex = d_vertices[index+1] - 1; // min value is stored in loc of last neighbour

		unsigned int val = d_segmentedMinScanOutput[minEdgeIndex];
		//unsigned int minWeight = val >> 22;
		unsigned int minVertex = val & (unsigned int)(pow(2.0,22)-1);

		d_successorArray[index] = minVertex;
	}
}

/*removes cycles from successor array*/
__global__ void eliminateCycles(unsigned int *d_successor, unsigned int noOfVertices) 
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfVertices) {
		unsigned int succIndex = d_successor[d_successor[index]];
		if(index == succIndex) {
			if(index < d_successor[index]) {
				d_successor[index] = index;
			 } else {
				d_successor[d_successor[index]]= d_successor[index];
			 }
		}
	}
}

/* hybrid implementation of markSegment function */
__global__ void markSegment1(unsigned int *d_edgeFlagArray, unsigned int *d_vertex, unsigned int noOfVertices)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

        if(index < noOfVertices && index > 0) {
                d_edgeFlagArray[d_vertex[index]] = 1;
        }
}

/*This function is to determine which edges are actually needed*/
__global__ void populatePArray(int *d_pickArray, unsigned int *d_vertices, unsigned int *d_successor, unsigned int *d_preIDs, unsigned int noOfVertices, unsigned int noOfEdges) 
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfEdges) {
		if(d_preIDs[index] != d_successor[d_preIDs[index]]) {
			if(d_preIDs[index] < (noOfVertices - 1))
				d_pickArray[index] = d_vertices[d_preIDs[index]+1] - 1;
			else
				d_pickArray[index] = noOfEdges - 1;	
		}
		else							
			d_pickArray[index] = -1;	
	}
}

/*This function determines which edges will be part of output*/
__global__ void AppendOutputEdges(int *d_pickArray, unsigned int * d_segmentedMinScanInput, unsigned int *d_segmentedMinScanOutput, unsigned int *d_MSTOutput, unsigned int *d_edgeMap, unsigned int noOfEdges)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfEdges && d_pickArray[index] >= 0) {
		unsigned int edgeid = d_edgeMap[index];
		unsigned int prev = 0;
		int temp = -1;
		unsigned int segmentedOutput = d_segmentedMinScanOutput[d_pickArray[index]];
		unsigned int currIndex = d_segmentedMinScanOutput[index];

		if(index > 0) {
			temp = d_pickArray[index-1];
			prev = d_segmentedMinScanOutput[index-1];	
		}

		if(d_pickArray[index] != temp) {
			if(currIndex == segmentedOutput) {
				d_MSTOutput[edgeid]=1;
			}		
		} else {
			if(currIndex != prev && currIndex == segmentedOutput) {	
				d_MSTOutput[edgeid]=1;
			}

		}
	}
}

/*This function sets each value of array equal to its index*/
__global__ void setIndices(unsigned int *d_arr,unsigned int size)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	if(index < size)
		d_arr[index] = index;
}

/* This function copies data from original successorArray so that it can be used for new computation*/
__global__ void setIndices1(unsigned int *d_arr,unsigned int size, int l)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
        if (index < size && index >= l)
                d_arr[index] = index;
}

__global__ void makeTempSuccCopy(unsigned int *d_successorArray, unsigned int* d_vertex, unsigned int *d_successorArrayTemp, unsigned int noOfVertices)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfVertices) {
		unsigned int t = d_successorArray[index];
		d_successorArrayTemp[index] = t;
	}
}

/* This function copies data from temporary successorArray so that it can be updated with correct value */
__global__ void updateSuccArray(unsigned int *d_successorArray, unsigned int* d_vertex, unsigned int *d_successorArrayTemp, unsigned int noOfVertices)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfVertices) {
		unsigned int t =  d_successorArrayTemp[index];
		d_successorArray[index] = t;
	}
}

/* This function uses pointer doubling to assign representative id to each vertex*/
__global__ void propagateRepVertexID(unsigned int *d_successorArray, bool *d_isSuccUpdated, unsigned int *d_previousIDs, unsigned int *d_successorArrayTemp, unsigned int noOfVertices)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfVertices) {
		unsigned int successor = d_successorArray[index];
		if(successor != d_successorArray[successor]) { //Eindex = 2 and end = 6 and  u = 2 and succ[u] = 2
			*d_isSuccUpdated=true;
			d_successorArrayTemp[index] = d_successorArray[successor];
		}
	}
}

/* This function iteratively sets s(s(u)) = u and propogates representative vertex id*/
void propagateID(unsigned int noOfBlocks_vertices, unsigned int noOfThreads_vertices)
{
	bool succchange;
	bool *d_isSuccUpdated;
        cudaMalloc(&d_successorArrayTemp, sizeof(int)*noOfVertices);
	cudaMalloc((void**)&d_isSuccUpdated, sizeof(bool));

	do
	{
		succchange=false;
		cudaMemcpy(d_isSuccUpdated, &succchange, sizeof(bool), cudaMemcpyHostToDevice);
		makeTempSuccCopy<<<noOfBlocks_vertices,noOfThreads_vertices>>>(d_successorArray, d_vertices, d_successorArrayTemp, noOfVertices);
		propagateRepVertexID<<<noOfBlocks_vertices,noOfThreads_vertices>>>(d_successorArray, d_isSuccUpdated, d_previousIDs,d_successorArrayTemp, noOfVertices);
		updateSuccArray<<<noOfBlocks_vertices,noOfThreads_vertices>>>(d_successorArray, d_vertices, d_successorArrayTemp, noOfVertices);
		cudaMemcpy(&succchange, d_isSuccUpdated, sizeof(bool), cudaMemcpyDeviceToHost);
	}while(succchange);

	cudaFree(d_successorArrayTemp);
	cudaFree(d_isSuccUpdated);
}

/*This function creates scan flag*/
void __global__ createScanFlag(unsigned int *d_vertexFlagArray, unsigned int *d_successorArray, unsigned int noOfVertices)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

        if(index < noOfVertices && index > 0)
	{
		unsigned int prev_val = d_successorArray[index-1];
		unsigned int curr_val = d_successorArray[index];
		if (prev_val != curr_val) {
			d_vertexFlagArray[index] = 1;
		}
	}	
}

/*This function assigns supervertex id to each vertex*/
__global__ void assignSuperVertexID(unsigned int *d_superVertex, unsigned int *d_indices, unsigned int *d_vertexFlagArray,unsigned int *d_previousIDs,unsigned int noOfVertices)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
	if(index < noOfVertices) {
		d_vertexFlagArray[d_indices[index]] = d_superVertex[index]; 
	}
}

/* This function updates supervertexid */
__global__ void updateSuperVertexID(unsigned int *d_superVertex,unsigned int *d_arr,unsigned int *d_vertexFlagArray, unsigned int noOfVertices)
{
        unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
        if(index < noOfVertices) {
        	unsigned int newId = d_vertexFlagArray[index];
                d_superVertex[index] = newId;
                
         }
}

/* This function removes self edges after successor array is created */
__global__ void removeSelfEdges(unsigned int *d_edges, unsigned int *d_prevIds,unsigned int *d_superVertexID, unsigned int noOfEdges)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index < noOfEdges) {
		 unsigned int uid = d_superVertexID[d_prevIds[index]]; //vause d_prevIds[index] is 1 to 6 but we need 0 to 5
		 unsigned int vid = d_superVertexID[d_edges[index]];
		 if(uid == vid) {
			d_edges[index]=INT_MAX; 
		 }
	}
}

/* This function is to assign new super edge id*/
__global__ void assignSuperEdgeId(unsigned int *d_superEdgeId, unsigned int *d_previousIds, unsigned int *d_superVertexId, unsigned int *d_edge, unsigned int noOfEdges)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

        if(index < noOfEdges)
        {
		unsigned int x = d_previousIds[index];
		unsigned int id = INT_MAX;
		if (x != INT_MAX && d_edge[index] != INT_MAX) {
			id = d_superVertexId[x];
		}
		d_superEdgeId[index] = id;		
	}
}


/* This function is to compress the edge list*/
__global__ void edgeCompression(unsigned int *d_edges, unsigned int *d_weights, unsigned int *d_vertex, unsigned int *d_segmentedMinScanInput, unsigned int *d_segmentedMinScanOutput, unsigned int *d_superVertexID, unsigned int *d_edgeMap, unsigned int *d_edgeMapCopy, unsigned int *d_edgeFlagArray, unsigned int *d_superEdgeId, unsigned int * d_edgeIndices, int *d_pickArray, unsigned int *d_size, unsigned int *d_edgeListSize, unsigned int *d_vertexListSize)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

        if(index < *d_size) {	
		unsigned int id = d_edgeIndices[index];
		if(d_superEdgeId[index] != INT_MAX && d_edges[id] != INT_MAX) {
			if(index == *d_size-1) {
				*d_edgeListSize = index + 1;
				*d_vertexListSize = d_superEdgeId[index] + 1;
			}

			d_segmentedMinScanOutput[index] = d_weights[id];
			d_segmentedMinScanInput[index] = d_superVertexID[d_edges[id]];
			d_pickArray[index] = d_superEdgeId[index];			
			d_edgeMapCopy[index] = d_edgeMap[id];
		}
	}

}

/*This function copies the temporary array to arrays which will be actually used*/
__global__ void copyArrays(unsigned int *d_edges, unsigned int *d_weights, unsigned int *vertex, unsigned int *d_segmentedMinScanInput, unsigned int *d_segmentedMinScanOutput, unsigned int *d_edgeMap, unsigned int *d_edgeMapCopy, unsigned int *d_edgeCopy, unsigned int *d_size)
{
	unsigned int index = blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

        if(index < *d_size) { 
        
        	unsigned int p = d_segmentedMinScanInput[index]; 
		d_edges[index] = p;
		unsigned int wt = d_segmentedMinScanOutput[index];
		d_weights[index] = wt;
		unsigned int mapVal =  d_edgeMapCopy[index];
		d_edgeMap[index] = mapVal;
	}
	
}

/*This function determines the new edge list*/
__global__ void makeEdgeList(unsigned int *d_edgeFlagArray, unsigned int *d_edges, unsigned int *d_superEdgeId, unsigned int *d_size, unsigned int noOfEdges)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index == 0) {
		d_edgeFlagArray[index] = 1;
	} else if(index < noOfEdges && index > 0) {
		if(d_superEdgeId[index-1] != INT_MAX && d_superEdgeId[index] == INT_MAX) {
			*d_size = index;
		}
		if(d_superEdgeId[index] > d_superEdgeId[index-1]) {
			d_edgeFlagArray[index] = 1;
		}			
        }
}

/*This function helps in creating new vertices list for next iteration*/
__global__ void CreateVertexListFlag(unsigned int *d_edgeFlagArray, unsigned int *d_vertices, int *d_pickArray, unsigned int noOfEdges)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;

	if(index == 0) {
		d_edgeFlagArray[index] = 1;
	} else if(index < noOfEdges && index > 0) {
		if(d_pickArray[index] > d_pickArray[index-1]) {
                	d_edgeFlagArray[index] = 1;
		}       
        }	
}

/*This function helps to build new vertex list*/
__global__ void BuildVertexList(unsigned int *d_vertices, unsigned int *d_edges, int *d_pickArray, unsigned int *d_edgeFlagArray, unsigned int noOfEdges)
{
	unsigned int index= blockIdx.x * NO_OF_THREADS_PER_BLOCK + threadIdx.x;
        if(index < noOfEdges && d_edgeFlagArray[index] == 1) {
		d_vertices[d_pickArray[index]] = index;
	}
}

/* Parse the input file to setup our graph
 * we set the relevant arrays here 
 */
void parseInputFile(char *fileName)
{
	unsigned int x,temp;
	unsigned int edgeNo, weightOfEdge;
	FILE *fp;
	fp = fopen(fileName,"r");
	
	printf("\n Parsing Input File: \n");
	fscanf(fp,"%d",&noOfVertices);
	vertices = (unsigned int *)malloc(sizeof(unsigned int) * noOfVertices);
	int i;
	for (i=0; i<noOfVertices; i++) {
		fscanf(fp,"%d %d",&x, &temp);
		vertices[i] = x;
	}

	fscanf(fp,"%d",&temp);
	fscanf(fp,"%d",&noOfEdges);
	
	edges = (unsigned int *)malloc(sizeof(unsigned int)*noOfEdges);
	weights = (unsigned int *)malloc(sizeof(unsigned int)*noOfEdges);

	for(i=0; i<noOfEdges; i++) {
		fscanf(fp,"%d %d",&edgeNo, &weightOfEdge);
		edges[i] = edgeNo;
		weights[i] = weightOfEdge;
	}
	
	printf("No. of Vertices in Input: %d\n",noOfVertices);
	printf("No. of Edges in Input: %d\n", noOfEdges);
	fclose(fp);
}

/* this is to setup configuration parameters for various primitives*/
void setupPlan()
{
	cudppCreate(&theCudpp);
	scan_config.algorithm = CUDPP_SCAN;
	scan_config.op = CUDPP_ADD;
	scan_config.datatype = CUDPP_UINT;
	scan_config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;

	segmented_min_scan_config.algorithm = CUDPP_SEGMENTED_SCAN;
	segmented_min_scan_config.op = CUDPP_MIN;
	segmented_min_scan_config.datatype = CUDPP_UINT;
	segmented_min_scan_config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;


	config_sort.algorithm = CUDPP_SORT_RADIX;
	config_sort.datatype = CUDPP_UINT;
	config_sort.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_KEY_VALUE_PAIRS;

	f = 0.05;
}

/* Dynamically allocate necessary arrays*/
void mallocArr()
{
 	cudaMalloc(&d_segmentedMinScanInput, sizeof(unsigned int )*noOfEdges);
        cudaMalloc(&d_weights, sizeof(unsigned int  )*noOfEdges);
        cudaMalloc(&d_edges, sizeof(unsigned int )*noOfEdges);
        cudaMalloc(&d_vertices, sizeof(unsigned int )*noOfVertices);
        cudaMalloc(&d_edgeFlagArray, sizeof(unsigned int )*noOfEdges);
        cudaMalloc(&d_segmentedMinScanOutput, sizeof(unsigned int )*noOfEdges);
        cudaMalloc(&d_successorArray, sizeof(unsigned int )*noOfVertices);
        cudaMalloc(&d_previousIDs, sizeof(unsigned int )*noOfEdges);
        cudaMalloc(&d_pickArray, sizeof(int )*noOfEdges);
	cudaMalloc(&d_superVertexID, sizeof(unsigned int )*noOfVertices);
        cudaMalloc(&d_MSTOutput, sizeof(unsigned int )*noOfEdges);
		
	cudaMalloc(&d_indices, sizeof(unsigned int )*noOfVertices);
	cudaMalloc(&d_vertexFlagArray, sizeof(unsigned int )*noOfVertices);
	cudaMalloc(&d_superVertexID, sizeof(unsigned int )*noOfVertices);
	cudaMalloc(&d_size, sizeof(unsigned int ));
 	cudaMalloc(&d_superEdgeId, sizeof(unsigned int )*noOfEdges);
	cudaMalloc(&d_edgeIndices, sizeof(unsigned int )*noOfEdges);
	cudaMalloc(&d_edgeListSize, sizeof(unsigned int ));
	cudaMalloc(&d_vertexListSize, sizeof(unsigned int ));
	cudaMalloc(&d_edgeMapCopy, sizeof(unsigned int )*noOfEdges);
	cudaMalloc(&d_edgeMap, sizeof(unsigned int )*noOfEdges);
	h_MSTOutput = (unsigned int *)malloc(sizeof(unsigned int )*noOfEdges);
}

/*Free the dynamically allocated memory. Do other cleanup here*/
void cleanUp()
{
	cudaFree(d_edgeIndices);
	cudaFree(d_superEdgeId);
	cudaFree(d_edgeMap);
	cudaFree(d_edgeMapCopy);
	cudaFree(d_superVertexID);		
	cudaFree(d_vertexFlagArray);
	cudaFree(d_indices);
	cudaFree(d_MSTOutput);
	cudaFree(d_previousIDs);
	cudaFree(d_pickArray);
	cudaFree(d_successorArray);
	cudaFree(d_segmentedMinScanOutput);
	cudaFree(d_edgeFlagArray);
	cudaFree(d_vertices);
	cudaFree(d_edges);
	cudaFree(d_weights);
	cudaFree(d_segmentedMinScanInput);
	cudaFree(d_size);
	cudaFree(d_edgeListSize);
	cudaFree(d_vertexListSize);
	
	cudppDestroy(theCudpp);
	free(h_MSTOutput);
	free(edges);
	free(vertices);
	free(weights);
}

/* Do basic initialization*/
void initialize()
{
	unsigned int i;

	cudaMemcpy(d_vertices, vertices, sizeof(unsigned int)*noOfVertices, cudaMemcpyHostToDevice);
	cudaMemcpy(d_edges, edges, sizeof(unsigned int)*noOfEdges, cudaMemcpyHostToDevice);
	cudaMemcpy(d_weights, weights, sizeof(unsigned int)*noOfEdges, cudaMemcpyHostToDevice);

	unsigned int *temp = (unsigned int *)malloc(sizeof(unsigned int)*noOfEdges);
	for(i=0; i<noOfEdges; i++)
		temp[i] = 0;
	cudaMemcpy(d_MSTOutput, temp, sizeof(unsigned int )*noOfEdges, cudaMemcpyHostToDevice);
	
	for(i=0; i<noOfEdges; i++)
		temp[i]=i;
	cudaMemcpy(d_edgeMap, temp, sizeof(unsigned int)*noOfEdges, cudaMemcpyHostToDevice);

	free(temp);
}

/* Helper function to determine no of threads to be used */
unsigned int getNoOfThreads(unsigned int size) {
	unsigned int threadsPerBlock;
	if (size <= 1024)
		threadsPerBlock = size;
	else
		threadsPerBlock = 1024;
	return threadsPerBlock;
}


void boruvka()
{
	int t;
	
	unsigned int noOfThreads_edge = getNoOfThreads(noOfEdges);
	unsigned int noOfBlocks_edge = (noOfEdges+1024)/noOfThreads_edge;
	unsigned int noOfThreads_vertices = getNoOfThreads(noOfVertices);
	unsigned int noOfBlocks_vertices = (noOfVertices+1024)/noOfThreads_vertices;
	cudaError_t error;	
	
	mergeEdgeAndWeight<<<noOfBlocks_edge, noOfThreads_edge>>>(d_segmentedMinScanInput, d_vertices, d_weights, d_edges, noOfEdges);

	error = cudaGetLastError();
       	if(error != cudaSuccess)
        {
                printf("0.1 CUDA error: %s\n", cudaGetErrorString(error));
                exit(-1);
        }

	t = noOfEdges * f;
 	if (noOfEdges >= 200)
        {
                unsigned int *temp_h_efa = (unsigned int *)malloc(t*sizeof(unsigned int));
		initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, noOfEdges,t);
               	int i;
		#pragma omp parallel for	
                for(i = 0; i<t;i++)
                        temp_h_efa[i] = 0;
		#pragma omp barrier
                cudaThreadSynchronize();
               	cudaMemcpy(d_edgeFlagArray, temp_h_efa, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                free(temp_h_efa);
        }
        else
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, noOfEdges);

	cudaThreadSynchronize();
	error = cudaGetLastError();
       	if(error != cudaSuccess)
        {
                printf("At line 613 CUDA error: %s\n", cudaGetErrorString(error));
                exit(-1);
        }

	markSegment<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_edgeFlagArray, d_vertices, d_edges, noOfVertices);
	
        error = cudaGetLastError();
        if(error != cudaSuccess)
       	{
        	printf("3  CUDA error: %s\n", cudaGetErrorString(error));
          	exit(-1);
       	}

	cudppPlan(theCudpp, &segmentedScanPlan_min,segmented_min_scan_config, noOfEdges, 1, 0 ); //Make the segmented min scan plan
	cudppSegmentedScan(segmentedScanPlan_min, d_segmentedMinScanOutput, d_segmentedMinScanInput, (const unsigned int *)d_edgeFlagArray, noOfEdges);
	cudppDestroyPlan(segmentedScanPlan_min);
		
	 error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("CUDA error: %s\n", cudaGetErrorString(error));
               // exit(-1);
        }
	
	createSuccArray<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_successorArray, d_vertices, d_segmentedMinScanInput, d_segmentedMinScanOutput, noOfVertices, noOfEdges);
	
	eliminateCycles<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_successorArray, noOfVertices);

 	t = noOfEdges * f;
        if (noOfEdges >= 200)
        {
               	unsigned int *temp_h_efa = (unsigned int *)malloc(t*sizeof(unsigned int));
                initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, noOfEdges,t);
                int i;
		#pragma omp parallel for
               	for(i = 0; i<t;i++)
                        temp_h_efa[i] = 0;
		#pragma omp barrier
                cudaThreadSynchronize();
                cudaMemcpy(d_edgeFlagArray, temp_h_efa, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
               	free(temp_h_efa);
        }
	else
            	initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, noOfEdges);

	cudaThreadSynchronize();
	markSegment1<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_edgeFlagArray, d_vertices, noOfVertices);
	
	cudppPlan(theCudpp, &scanPlan, scan_config, noOfEdges, 1, 0);
	cudppScan(scanPlan, d_previousIDs, d_edgeFlagArray, noOfEdges);
	cudppDestroyPlan(scanPlan);
	
	 error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("At line 668 CUDA error: %s\n", cudaGetErrorString(error));
        }

 	t = noOfEdges * f;
        if(noOfEdges >= 200)
        {
                unsigned int *temp_h_pa = (unsigned int *)malloc(t*sizeof(unsigned int));
        	initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>((unsigned int*)d_pickArray, noOfEdges, t);

                int i;
		#pragma omp parallel for
                for(i = 0; i<t;i++)
                        temp_h_pa[i] = 0;
		#pragma omp barrier
                cudaThreadSynchronize();
                cudaMemcpy(d_pickArray, temp_h_pa, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                free(temp_h_pa);
        }
	else
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>((unsigned int*)d_pickArray, noOfEdges);

	cudaThreadSynchronize();

	populatePArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_pickArray, d_vertices, d_successorArray, d_previousIDs, noOfVertices, noOfEdges);
	
	AppendOutputEdges<<<noOfBlocks_edge, noOfThreads_edge>>>(d_pickArray, d_segmentedMinScanInput, d_segmentedMinScanOutput, d_MSTOutput, d_edgeMap, noOfEdges);

	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("At line 698  CUDA error: %s\n", cudaGetErrorString(error));
        }
	
	propagateID(noOfBlocks_vertices, noOfThreads_vertices);
      
        t = noOfVertices*f;

	if(noOfVertices >= 20)
	{
		unsigned int *temp_h_setIndices = (unsigned int *)malloc(t*sizeof(unsigned int));
	        setIndices1<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_indices, noOfVertices, t);
		int i;
		#pragma omp parallel for
		for(i = 0; i<t;i++)
			temp_h_setIndices[i] = i;
		#pragma omp barrier
		cudaThreadSynchronize();
		cudaMemcpy(d_indices, temp_h_setIndices, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
		free(temp_h_setIndices);		
	}
	else
        	setIndices1<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_indices, noOfVertices, 0);

        cudaThreadSynchronize();

	cudppPlan(theCudpp, &sortPlan, config_sort, noOfVertices, 1, 0);
        cudppRadixSort(sortPlan, d_successorArray, d_indices, noOfVertices);
        cudppDestroyPlan(sortPlan);

	initArray<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_vertexFlagArray,noOfVertices);
	createScanFlag<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_vertexFlagArray, d_successorArray,noOfVertices);
	
	cudppPlan(theCudpp, &scanPlan, scan_config, noOfVertices, 1, 0);
        cudppScan(scanPlan, d_superVertexID, d_vertexFlagArray, noOfVertices);
        cudppDestroyPlan(scanPlan);
	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("At line 736 CUDA error: %s\n", cudaGetErrorString(error));
        }	
	
	assignSuperVertexID<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_superVertexID,d_indices,d_vertexFlagArray,d_previousIDs,noOfVertices);
       
	updateSuperVertexID<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_superVertexID,d_indices,d_vertexFlagArray, noOfVertices);
	removeSelfEdges<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edges,d_previousIDs,d_superVertexID,noOfEdges);
	
	assignSuperEdgeId<<<noOfBlocks_edge, noOfThreads_edge>>>(d_superEdgeId,d_previousIDs, d_superVertexID, d_edges, noOfEdges);
 
 	t = noOfEdges*f;
        //printf("noOfVertices = %d and point  = %d\n",noOfVertices, t);
	if (noOfEdges >= 200)
        {
                unsigned int *temp_h_setIndices = (unsigned int *)malloc(t*sizeof(unsigned int));
                setIndices1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeIndices, noOfEdges, t);
               	int i;
		#pragma omp parallel for
                for(i = 0; i<t;i++)
                        temp_h_setIndices[i] = i;
		#pragma omp barrier
                cudaThreadSynchronize();
               	cudaMemcpy(d_edgeIndices, temp_h_setIndices, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
        	free(temp_h_setIndices);
	}
       	else
                setIndices1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeIndices,noOfEdges,0);

	cudaThreadSynchronize();
	cudppPlan(theCudpp, &sortPlan, config_sort, noOfEdges, 1, 0);
        cudppRadixSort(sortPlan, d_superEdgeId, d_edgeIndices, noOfEdges);
        cudppDestroyPlan(sortPlan);

 	t = noOfEdges * f;
        if (noOfEdges >= 200)
        {
                unsigned int *temp_h_efa = (unsigned int *)malloc(t*sizeof(unsigned int));
                initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>((unsigned int*)d_edgeFlagArray, noOfEdges, t);
                int i;
		#pragma omp parallel for
                for(i = 0; i<t;i++)
                       	temp_h_efa[i] = 0;
		#pragma omp barrier
                cudaThreadSynchronize();
                cudaMemcpy(d_edgeFlagArray, temp_h_efa, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                free(temp_h_efa);
        }
        else
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray,noOfEdges);

	cudaThreadSynchronize();
	unsigned int  h_size = noOfEdges + 1;
	cudaMemcpy(d_size,&h_size,sizeof(unsigned int ), cudaMemcpyHostToDevice); 
	makeEdgeList<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, d_edges, d_superEdgeId, d_size, noOfEdges);	
	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("At line 793  CUDA error: %s\n", cudaGetErrorString(error));
        }

	unsigned int  zero = 0;
	cudaMemcpy(d_edgeListSize, &zero, sizeof(unsigned int ), cudaMemcpyHostToDevice);
	cudaMemcpy(d_vertexListSize, &zero, sizeof(unsigned int ), cudaMemcpyHostToDevice);

	t = noOfEdges * f;
       	if (noOfEdges >= 200)
        {
                unsigned int *temp_arr = (unsigned int *)malloc(t*sizeof(unsigned int));
                initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_segmentedMinScanInput, noOfEdges,t);
		initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_segmentedMinScanOutput, noOfEdges, t);
                initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>((unsigned int*)d_pickArray, noOfEdges, t);
                initArray1<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeMapCopy, noOfEdges, t);      

         	int i;
		#pragma omp parallel for
               	for(i = 0; i<t;i++)
                       	temp_arr[i] = 0;
		#pragma omp barrier
                cudaThreadSynchronize();
                cudaMemcpy(d_segmentedMinScanInput, temp_arr, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                cudaMemcpy(d_segmentedMinScanOutput, temp_arr, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                cudaMemcpy(d_pickArray, temp_arr, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);
                cudaMemcpy(d_edgeMapCopy, temp_arr, sizeof(unsigned int )*t, cudaMemcpyHostToDevice);

		free(temp_arr);
        }
        else
	{
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_segmentedMinScanInput, noOfEdges);
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_segmentedMinScanOutput, noOfEdges);
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>((unsigned int*)d_pickArray, noOfEdges);
		initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeMapCopy, noOfEdges);
	
	}
	
	cudaThreadSynchronize();
	cudaMemcpy(&h_size,d_size,sizeof(unsigned int ), cudaMemcpyDeviceToHost);

	unsigned int  noOfThreads_new = getNoOfThreads(h_size);
        unsigned int  noOfBlocks_new = (h_size+1024)/noOfThreads_new;

	edgeCompression<<<noOfBlocks_new, noOfThreads_new>>>(d_edges, d_weights, d_vertices, d_segmentedMinScanInput, d_segmentedMinScanOutput, d_superVertexID, d_edgeMap, d_edgeMapCopy, d_edgeFlagArray, d_superEdgeId, d_edgeIndices, d_pickArray, d_size, d_edgeListSize, d_vertexListSize);
	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("At line 841  CUDA error: %s\n", cudaGetErrorString(error));
        }
	copyArrays<<<noOfBlocks_new, noOfThreads_new>>>(d_edges, d_weights, d_vertices, d_segmentedMinScanInput, d_segmentedMinScanOutput, d_edgeMap, d_edgeMapCopy, d_edgeFlagArray, d_size);
	
	initArray<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, noOfEdges);
        initArray<<<noOfBlocks_vertices, noOfThreads_vertices>>>(d_vertices, noOfVertices);
	CreateVertexListFlag<<<noOfBlocks_edge, noOfThreads_edge>>>(d_edgeFlagArray, d_vertices, d_pickArray, noOfEdges);
	BuildVertexList<<<noOfBlocks_edge, noOfThreads_edge>>>(d_vertices, d_edges, d_pickArray, d_edgeFlagArray, noOfEdges);
	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("after build vertex listlast  CUDA error: %s\n", cudaGetErrorString(error));
        }

	cudaMemcpy(&noOfEdges, d_edgeListSize, sizeof(unsigned int ), cudaMemcpyDeviceToHost);
	cudaMemcpy(&noOfVertices, d_vertexListSize, sizeof(unsigned int ), cudaMemcpyDeviceToHost);

	printf("for next round, no of edges = %d and no of vertices = %d\n",noOfEdges, noOfVertices);	
 	error = cudaGetLastError();
        if(error != cudaSuccess)
        {
                printf("last  CUDA error: %s\n", cudaGetErrorString(error));
        }
}


int  main (int  argc, char** argv)
{
	unsigned int noOfMSTEdges = 0;
	unsigned long long int finalMSTWeight = 0;
	unsigned int i;
	parseInputFile(argv[1]);
	noOfVerticesOriginal = noOfVertices;
	noOfEdgesOriginal = noOfEdges;
	omp_set_dynamic(0);
        omp_set_num_threads(omp_get_num_procs());
	
	mallocArr();
	initialize();
	setupPlan();

	struct timeval  tv1, tv2;
        gettimeofday(&tv1, NULL);
	do {
		boruvka();
	}while(noOfVertices > 1);
	cudaThreadSynchronize();

	gettimeofday(&tv2, NULL);
        printf ("Total Execution time = %f seconds\n", (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
	cudaMemcpy(h_MSTOutput, d_MSTOutput, sizeof(unsigned int )*noOfEdgesOriginal, cudaMemcpyDeviceToHost);

	for(i=0; i<noOfEdgesOriginal; i++) {
		if(h_MSTOutput[i] == 1) {
			//printf("%d %d\n", edges[i], weights[i]);
			finalMSTWeight += weights[i];
			noOfMSTEdges++;
		}
	}
	printf("\nNo. of edges in MST [must be equal to (%d-1)]: %d\n", noOfVerticesOriginal, noOfMSTEdges);
	printf("Final Weight of resultant MST: %llu\n", finalMSTWeight);

	cleanUp();
	return 0;	
}

