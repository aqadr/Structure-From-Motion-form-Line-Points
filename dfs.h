#ifndef DFS_H
#define DFS_H



#include <iostream>
#include <list>


using namespace std;

class Graph
{
    int numV;   //number of vertices
    list<int> * adj; //a list of pointers. Where the pointers point to an arrray containing adjacency
    void DFSUtil(int v, bool visited[]);
public:
    Graph(int V);
    void addEdge(int v, int w);
    void DFS(int v);
};


//initializing the graph
Graph::Graph(int V)
{
    this->V=V;
    adj=new list<int>[V];
}

//add edge for each node in the graph
Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
}

void Graph::DFSUtil(int v, bool visited[])
{
    visited[v]=true;
    list<int>::iterator it;
    for(it=adj[v].begin();it!=adj[v].end();++it)
    {
        if(!visited[*it])
        {
            DFSUtil(*it,visited);
        }
    }
}

void Graph::DFS(int v)
{
    bool *visited = new bool[numV];
    for(int i=0;i<numV;++i)
    {
        visited[i]=false;
    }

    DFSUtil(v,visited);

}

#endif // DFS_H

