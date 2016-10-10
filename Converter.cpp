// Converter.cpp : Converts a DIMACS Graph to My Input Format
// We've tested it only for US Network graphs.
// Does not remove selfloops, negative weigts and does not make the graph connected.

#include <fstream>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;

struct node
{
	vector<int>edges;
	vector<int>weights;
};

node *Nodes;

int main(int argc, char* argv[])
{
	int no_of_nodes, no_of_edges;
	char ch;
	int x,y,weight;
	char str[250];
	int p=0;
	fstream f;
	f.open(argv[1],ios::in);
	if(!f)cout<<"Error reading file\n";
	while(ch!='a')
	{
	f>>ch;
	if(ch=='c')f.getline(str,250,'\n');
	if(ch=='p')
		{
			f>>ch>>ch;
			f>>no_of_nodes;
			f>>no_of_edges;
		}
	}

	Nodes = new node[no_of_nodes];

	f.seekg(-1,ios::cur);

	p=0;
	int selfloops=0,negativeweight=0;
	for(int i=0;i<no_of_edges;i++)
	{

		f>>ch>>x>>y>>weight;
		x--;
		y--;
		if(weight<0)
			{
			weight=-weight;
			negativeweight++;
			}
		if(x!=y)
		{
			Nodes[x].edges.push_back(y);
			Nodes[x].weights.push_back(weight);
			p++;
		}
		else
			selfloops++;

	}
	f.close();

	cout<<"Finished Reading file\n";
	cout<<"No_of_nodes"<<" "<<no_of_nodes<<" Edges	"<<no_of_edges<<" Total edges Added "<<p<<" Self Loops Not added "<<selfloops<<" Negative Weights "<<negativeweight<<endl; 

	bool found=false;

	for(int i=0;i<no_of_nodes;i++)
	{
		for(int j=0; j<Nodes[i].edges.size();j++)
		{
			bool notfound=false;
			int p = Nodes[i].edges[j];
			for(int k=0; k<Nodes[p].edges.size(); k++)
			{
				if( i == Nodes[p].edges[k])
					found=true;
			}
			if(!found)
				printf("No opp edge found for %d--->%d\n", i,p);
		}
	}

	cout<<"Finished checking Edge to Edge mask\n";

	fstream h;
	h.open(argv[2],ios::out);

	int starting=0;
	no_of_edges=0;

	h<<no_of_nodes<<"\n";
	
	for(int i=0;i<no_of_nodes;i++)
	{
	h<<starting<<" "<<Nodes[i].edges.size()<<"\n";
	starting+=Nodes[i].edges.size();
	no_of_edges+=Nodes[i].edges.size();
	}

	h<<"0\n";
	h<<no_of_edges<<"\n";

	for(int i=0;i<no_of_nodes;i++)
	for(int j=0;j<Nodes[i].edges.size();j++)
		h<<Nodes[i].edges[j]<<" "<<Nodes[i].weights[j]<<" ";

	h.close();

	delete []Nodes;
	return 0;
}

