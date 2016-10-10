#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include<stdlib.h>

using namespace std;

/*
 * A Generic and Highly Efficient Parallel Variant of Boruvka's Algorithm
 * http://ieeexplore.ieee.org/xpls/icp.jsp?arnumber=7092783&tag=1#search-basic
 */

void findMinEdgePerVertex(int *vertices, int *edges, int *weights, int *vertexMinEdge, int noOfVertices, int noOfEdges)
{
	int i, min_weight, min_edge, index;
	int j= 0;
	// for (i = 0; i < noOfVertices; i++) {
	// 	printf("iter vertices: %d - %d\n", i, vertices[i]);
	// }

	// for (i = 0; i < noOfEdges; i++) {
	// 	printf("iter edges: %d - %d\n", i, edges[i]);
	// }

	for (i = 0; i < noOfVertices; i++) {
		//printf("i = %d\n",i);
		j = vertices[i];
		// index = vertices[i];
		min_weight = weights[j];
		min_edge = edges[j];
		if (i == 4095)
			printf("j = %d\n", j);
		// if (i == 7) {
		// 	printf("%d - %d - %d\n", index, min_weight, min_edge);
		// }
		while (j!=noOfEdges && j != noOfEdges - 1 && (j < vertices[i + 1] || i + 1 == noOfVertices)) {
			/*if (i == 4095) {
				cout<<"j = "<<j<<endl;
				cout<<"noOfEdges = "<<noOfEdges<<endl;
				cout<<"diff = "<<noOfEdges - j<<endl;
				break;
				
			}*/
			if (weights[j] < min_weight) {
				min_weight = weights[j];
				min_edge = edges[j];
			}
			j++;
		}
		vertexMinEdge[i] = min_edge;
	}
}

void removeMirroredEdges(int *vertices, int *edges, int *weights, int *vertexMinEdge, int vertexId, int noOfVertices, int noOfEdges) {
	int i;

	if (vertexId == vertexMinEdge[vertexMinEdge[vertexId]]) {
		vertexMinEdge[vertexId] = -1;
	}
}

void insert_new_edges(int *vertices, int *edges, int *weights, int *color, int *successor, int *newVertices, int *newEdges, int *newWeights, int *newEdgeSize, int noOfVertices, int noOfEdges)
{
	int i, j = 0, k = 0, l = 0;
	int count, iterator = 1, vertexIter = 0;
	int startEdgeIndex, endEdgeIndex;

	while (vertexIter < noOfVertices) {
		count = 0;
		while (successor[vertexIter] != 1) {
			vertexIter++;
		}
		if (vertexIter >= noOfVertices) {
			break;
		}

		for (i = 0; i < noOfVertices; i++) {
			if (color[vertexIter] == color[i]) {
				startEdgeIndex = vertices[i];
				endEdgeIndex = i + 1 < noOfVertices ? vertices[i + 1] - 1 : noOfEdges - 1;
				if (newVertices[l] == -1) {
					newVertices[l] = k;
				}
				for (j = startEdgeIndex; j <= endEdgeIndex; j++) {
					if (color[vertexIter] != color[edges[j]]) {
						newEdges[k] = edges[j];
						newWeights[k] = weights[j];
						k++;
					}
				}
			}
		}
		vertexIter++;
		l++;
	}
	*newEdgeSize = k > 0? k - 1 : 0;
}

void parseGraph(int **vertices, int **edges, int **weights, int *noOfVertices, int *noOfEdges)
{
	int x, value, junk, i;
	fstream f;
	f.open("ipfile", ios::in);
	f >> x;
	*noOfVertices = x;
	*vertices = (int *)malloc(x * sizeof(int));
	for (i = 0; i < x; i++) {
		f >> value >> junk;
		//cout<<"value = "<<value<<endl;
		(*vertices)[i] = value;
	}

	f >> junk;
	f >> x;
	*noOfEdges = x;

	*edges = (int *)malloc(x * sizeof(int));
	*weights = (int *)malloc(x * sizeof(int));

	for (i=0; i<x; i++) {
		f >> value >>junk;
		//cout<<"junk = "<<junk<<endl;
		(*edges)[i] = value;
		(*weights)[i] = junk;
	}

}

int main(int argc, char *argv[])
{
	// int vertices[] = {0, 3, 6, 9, 12, 18, 21};
	// int edges[] = {1, 3, 4, 0, 2, 4, 1, 4, 6, 0, 4, 5, 0, 1, 2, 3, 5, 6, 3, 4, 6, 2, 4, 5};
	// int weights[] = {4, 12, 9, 4, 8, 7, 8, 5, 2, 12, 3, 6, 9, 7, 5, 3, 1, 13, 6, 1, 11, 2, 13, 11};
	// int vertices[] = {0, 2, 4, 7, 9, 11};
	// int edges[] = {1, 4, 0, 5, 3, 4, 5, 2, 5, 0, 2, 1, 2, 3};
	// int weights[] = {40, 9, 40, 11, 4, 10, 3, 4, 56, 9, 10, 11, 3, 56};
	// princeton 8-16

	int *vertices, *edges, *weights, noOfVertices, noOfEdges, i;

	parseGraph(&vertices, &edges, &weights, &noOfVertices, &noOfEdges);
	printf("%d\n", noOfEdges);

	// for (i = 0; i < noOfVertices; i++) {
	// 	printf("%d - %d\n", i, vertices[i]);
	// }

	// printf("*****\n");
	// for (i = 0; i < noOfEdges; i++) {
	// 	printf("%d - %d\n", edges[i], weights[i]);
	// }

	// return 0;



	// int vertices[] = {0, 4, 8, 13, 16, 20, 23, 27};
	// int edges[] = {2, 4, 6, 7, 2, 3, 5, 7, 0, 1, 3, 6, 7, 1, 2, 6, 0, 5, 6, 7, 1, 4, 7, 0, 2, 3, 4, 0, 1, 2, 4, 5};
	// int weights[] = {26, 38, 58, 16, 36, 29, 32, 19, 26, 36, 17, 40, 34, 29, 17, 52, 38, 35, 93, 37, 32, 35, 28, 58, 40, 52, 93, 16, 19, 34, 37, 28};
	// int noOfVertices = sizeof(vertices) / sizeof(vertices[0]);
	// int noOfEdges = sizeof(edges) / sizeof(edges[0]);
	
	int countEdges = 0;

	// int i,
	
	int colorIndex = 0;
	int vertexMinEdge[noOfVertices];
	int color[noOfVertices];
	int successor[noOfVertices];
	int flag[noOfVertices];
	int exclusivePrefixSum[noOfVertices];
	int colorCount;

	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	for (i = 0; i < noOfVertices; i++) {
		color[i] = -1;
		successor[i] = 0;
	}

	 printf("1\n");
	while (noOfVertices > 1) {
		//cout<<"1"<<endl;
		findMinEdgePerVertex(vertices, edges, weights, vertexMinEdge, noOfVertices, noOfEdges);
		 
		// printf("2\n");
		for (i = 0; i < noOfVertices; i++) {

			removeMirroredEdges(vertices, edges, weights, vertexMinEdge, i, noOfVertices, noOfEdges);
		}

		//printf("3\n");
		for (i = 0; i < noOfVertices; i++) {
			if (vertexMinEdge[i] != -1 && i < noOfVertices) {
				// printf("MST edge: %d - %d\n", i, vertexMinEdge[i]);
				countEdges++;
			}
		}

		colorCount = 0;
		for (i = 0; i < noOfVertices; i++) {
			if (vertexMinEdge[i] != -1) {
				if (color[i] == -1) {
					if (color[vertexMinEdge[i]] == -1) {
						successor[vertexMinEdge[i]] = 1;
						color[vertexMinEdge[i]] = colorIndex;
						colorIndex++;
						colorCount++;
					}
					color[i] = color[vertexMinEdge[i]];
				}
			}
		}

		exclusivePrefixSum[0] = 0;
		for (i = 1; i < noOfVertices; i++) {
			exclusivePrefixSum[i] = exclusivePrefixSum[i - 1] + successor[i - 1];
		}

		// printf("4\n");
		int *newVertices = (int *)malloc(colorCount*sizeof(int));
		int *newEdges;
		int *newWeights;

		int newVerticesSize = sizeof(newVertices) / sizeof(newVertices[0]);
		// printf("Total number of new vertices in this MST = %d\n", newVerticesSize);
		int newEdgeSize;

		// printf("5\n");
		if (newVerticesSize <= 1) {
			break;
		}

		for (i = 0; i < newVerticesSize; i++) {
			newVertices[i] = -1;
		}


		insert_new_edges(vertices, edges, weights, color, successor, newVertices, newEdges, newWeights, &newEdgeSize, noOfVertices, noOfEdges);

		// printf("6\n");
		// for (i = 0; i < noOfEdges; i++) {
		// 	printf("newEdges: %d - %d\n", i, edges[i]);
		// }

		newEdges = (int *)malloc(newEdgeSize* sizeof(int));
		newWeights = (int *)malloc(newEdgeSize * sizeof(int));
		// printf("noOfEdges - %d\n", newEdgeSize);
		memcpy(vertices, newVertices, newVerticesSize* sizeof(int));
		// printf("61\n");
		// printf("%d\n", );
		memcpy(edges, newEdges, newEdgeSize* sizeof(int));
		// printf("62\n");
		memcpy(weights, newWeights, newEdgeSize* sizeof(int));
		// printf("7\n");
		noOfVertices = newVerticesSize;
		noOfEdges = newEdgeSize;
		// printf("8\n");
		free(newVertices);
		free(newEdges);
		free(newWeights);
	}
	printf("Total number of edges in this MST = %d\n", countEdges);
	// printf("Total number of vertices in this MST = %d\n", noOfVertices);
	gettimeofday(&tv2, NULL);
	printf ("Total Execution time = %f seconds\n", (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
	free(vertices);
	free(edges);
	free(weights);
}



