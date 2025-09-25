#ifndef GRAPH_H
#define GRAPH_H
#include "Vertex.h"
#include "Edge.h"

#include <vector>
#include <string>

template <typename T>
class Graph {
public:
    Graph() {}

    void insert_vertex(const Vertex<T>& ver);
    void add_edge(const Vertex<T>& ver1, const Vertex<T>& ver2, int dist, int cost); //connect ver1 with ver2
    void print() const;

    void shortest_path(const Vertex<T>& src, const Vertex<T>& dest);
    void shortest_paths_to_state(const Vertex<T>& origin, const std::string& dest_state, const std::vector<std::pair<std::string, std::string>>& airport_to_state);
    void shortest_path_with_stops(const Vertex<T>& src, const Vertex<T>& dest, int n_stops);
    void count_direct_connections();
    Graph<T> create_undirected_graph();
    void prim_mst();
    void kruskal_mst(); 
private:
    std::vector<Vertex<T>> vertices; //nodes
    std::vector<std::vector<Edge>> edges; //connections
    int get_vertex_index(const Vertex<T>& ver);
    void clean_visited();
};

#endif
