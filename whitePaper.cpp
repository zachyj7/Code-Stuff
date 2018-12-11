#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

// A utility function to print the adjacency list representation of graph 
void printGraph(vector<vector<int> > adj) 
{ 
    for (int v = 0; v < adj.size(); ++v) 
    { 
        cout << "\nAdjacency list of vertex " << v << "\n head "; 
        for (int i = 0; i < adj[v].size(); ++i) 
		{
           cout << " " << adj[v][i]; 
		}
        cout << endl; 
    } 
} 

// Function that creates the Adjecency List from the input file
vector<vector<int> > createGraph(char*& graphFileName, int maxNode) {
	vector<int> innerAdj;
	vector<vector<int> > adj(maxNode + 1, innerAdj);
	ifstream source;
	string line;
	int Node1;
	int Node2;
    source.open(graphFileName);
	stringstream ss;  
	while (!source.eof())
    {
		getline(source, line);
		ss << line;
		ss >> Node1 >> Node2;
		adj[Node1].push_back(Node2);
    }
	return adj;
}

// Gets the total number of nodes in a directed graph file
int getNodeCount(char*& graphFileName) {
	string line;
	ifstream source;
	int maxNode = 0;
	int Node1;
	int Node2;
    source.open(graphFileName);
	stringstream ss;  
	while (getline(source, line))
    {
		ss << line;
		ss >> Node1 >> Node2;
		if (Node1 > maxNode) {
			maxNode = Node1;
		}
    }
	return maxNode;
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int K = atoi(argv[2]);
	float D = atof(argv[3]);
	int p = 4;
	
	// Variables used for calculating the top 5 highest pageranks
	int first = 0, second = 0, third = 0, fourth = 0, fifth = 0;
	int firstIndex = 0, secondIndex = 0, thirdIndex = 0, fourthIndex = 0, fifthIndex = 0;
	
	char* graphFileName(argv[1]);
	int maxNode = getNodeCount(graphFileName);
	
	vector<vector<int> > adj = createGraph(graphFileName, maxNode);
    //printGraph(adj); 
	vector<int> adjVisits(maxNode + 1, 0);	// Create the visit counter array

	//-------------------------------Start parallel section---------------------------//
	
	omp_set_num_threads(p);
	double dampingRatio = D * 100;	
	int rankTimes = (adj.size() / p);
	//cout << "Each rank does " << rankTimes << " nodes and the full array is this big: " << adj.size() <<endl;
	double time = omp_get_wtime();
	#pragma omp parallel
	{
		int randNum;
		int randNode;
		int currentNode;
		int listSize;
		int index;
		int rank = omp_get_thread_num();
		
		// Get the section start and finish for each rank
		int sectionStart = ((rank + 1) * rankTimes) - rankTimes;
		int sectionEnd = ((rank + 1) * rankTimes) - 1;
		
		for (int i = sectionStart; i <= sectionEnd; i++) {
			currentNode = i;
			for (int k = 0; k < K; k++) {
				listSize = adj[currentNode].size();
				randNum = rand() % 100 + 1;		// Flip coin
				if (randNum < dampingRatio) {	// If heads, choose random node
					randNode = rand() % (maxNode + 1);
					currentNode = randNode;
					adjVisits[randNode]++;		// Update the node visit count
				} else {	
					if (listSize > 0) {			// If tails, choose random adjacent member
						index = rand() % listSize;
						if (adj.at(currentNode).at(index) < adj.size()) {
							currentNode = adj[currentNode][index];
							adjVisits[currentNode]++;	// Update the node visit count
						}
					}
				}
			}
		}
		
	}
	time = omp_get_wtime() - time;
	
	//------------------------End parallel section----------------------------//
	
	
	// Figures out the top 5 highest page ranks
	for (int j = 0; j < adjVisits.size(); j++) {
		if (adjVisits[j] > first) {
			fifth = fourth;
			fourth = third;
			third = second;
			second = first;
			first = adjVisits[j];
			fifthIndex = fourthIndex;
			fourthIndex = thirdIndex;
			thirdIndex = secondIndex;
			secondIndex = firstIndex;
			firstIndex = j;
			
		} else if (adjVisits[j] > second) {
			fifth = fourth;
			fourth = third;
			third = second;
			second = adjVisits[j];
			fifthIndex = fourthIndex;
			fourthIndex = thirdIndex;
			thirdIndex = secondIndex;
			secondIndex = j;
			
		} else if (adjVisits[j] > third) {
			fifth = fourth;
			fourth = third;
			third = adjVisits[j];
			fifthIndex = fourthIndex;
			fourthIndex = thirdIndex;
			thirdIndex = j;
		} else if (adjVisits[j] > fourth) {
			fifth = fourth;
			fourth = adjVisits[j];
			
			fifthIndex = fourthIndex;
			fourthIndex = j;
		} else if (adjVisits[j] > fifth) {
			fifth = adjVisits[j];
			
			fifthIndex = j;
		}
		
		
	}
	
	cout << "Top 5 Nodes for dependencies in order are: " << endl;
	cout << "#1 Node Number " << firstIndex << " had " << first << " dependencies" << endl; 
	cout << "#2 Node Number " << secondIndex << " had " << second << " dependencies" << endl;
	cout << "#3 Node Number " << thirdIndex << " had " << third << " dependencies" << endl;
	cout << "#4 Node Number " << fourthIndex << " had " << fourth << " dependencies" << endl;
	cout << "#5 Node Number " << fifthIndex << " had " << fifth << " dependencies" << endl;
	cout << "Time taken for parallel region: " << time << " seconds" << endl;
	cout << endl;
	return 0;
}