#include "Graph.h"
#include "MinHeap.h"
#include "MinHeap.cpp"
#include "DisjointSet.h"

#include <string>
#include <unistd.h> 
#include <iostream>
#include <iomanip> 

#define INT_MAX 1000


template <typename T>
void Graph<T>::insert_vertex(const Vertex<T>& ver) {
    if (get_vertex_index(ver) == -1) {
        vertices.push_back(ver); //insert the vertex to the array of vertices
        std::vector<Edge> tmp;
        edges.push_back(tmp); //insert empty vector to the edges
    }
}

template <typename T>
int Graph<T>::get_vertex_index(const Vertex<T>& ver) {
    for(int i = 0; i < vertices.size(); i++) {
        if (vertices[i].getData() == ver.getData()) {
            return i;
        }
    }
    return -1;
}

template <typename T>
void Graph<T>::add_edge(const Vertex<T>& ver1, const Vertex<T>& ver2, int dist, int cost) {
    int i1 = get_vertex_index(ver1);
    int i2 = get_vertex_index(ver2);
    if (i1 == -1 || i2 == -1) {
        throw std::string("Add_edge: incorrect vertices");
    }
    Edge v(i1, i2, dist, cost);
    Edge v2(i2, i1, dist, cost);
    edges[i1].push_back(v);
    if (i1 != i2) {
        edges[i2].push_back(v2);
    }
}

template <typename T>
void Graph<T>::print() const {
    for (int i = 0; i < vertices.size(); i++) {
        std::cout << "{ " << vertices[i].getData() << ": ";
        for(int j = 0; j < edges[i].size(); j++) {
            std::cout << '{' << vertices[edges[i][j].dest].getData() << ", ";
            std::cout << "Distance: " << edges[i][j].distance << ", Cost: " << edges[i][j].cost << "} ";
        }
        std::cout << " }\n";
    }
}

template <typename T>
void Graph<T>::clean_visited() {
    for(Vertex<T>& v : vertices) {
        v.setVisited(false);
    }
}

//Task 2: Find the shortest path between the given origin airport and destination airport. Dijkstra's algorithm
template<typename T>
void Graph<T>::shortest_path(const Vertex<T>& src, const Vertex<T>& dest) {
    int i_src = get_vertex_index(src);
    int i_dest = get_vertex_index(dest);

    //If the origin or the destination airport does not exists, return
    if (i_src == -1 || i_dest == -1) {
        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData() << ": None\n";
        return;
    }

    //Initialize all edges to unvisited
    clean_visited();

    //Initialize
    std::vector<int> distances(vertices.size(), INT_MAX);
    std::vector<int> costs(vertices.size(), INT_MAX);
    std::vector<int> predecessors(vertices.size(), -1);
    distances[i_src] = 0;
    costs[i_src] = 0;

    //Start from the source vertix
    MinHeap<Edge> heap;
    int vertices_visited = 0;
    int cur_ver = i_src;

    //While have not visit all vertices
    while (vertices_visited < vertices.size()) {
        int i = cur_ver;
        //Go through each neighbor of the current vertext that is unvisited
        for (int j = 0; j < edges[i].size(); j++) {
            int i_adjacent_ver = edges[i][j].dest;

            if (!vertices[i_adjacent_ver].getVisited()) {
                heap.insert(edges[i][j]);
                int new_distance = distances[i] + edges[i][j].distance;
                
                //If the distance of a neighbor is shorter than the previous neighbor, update the path
                if (new_distance < distances[i_adjacent_ver]) {
                    distances[i_adjacent_ver] = new_distance;
                    costs[i_adjacent_ver] = costs[i] + edges[i][j].cost;
                    predecessors[i_adjacent_ver] = i;

                    
                }
            }
        }

        //The shortest path will be store until the shortest found
        if (heap.is_empty()) {
            break;
        }
        Edge e = heap.delete_min();
        cur_ver = e.dest;
        vertices[cur_ver].setVisited(true);
        vertices_visited++;
    }

    clean_visited();
    // Output the path
    if (distances[i_dest] == INT_MAX) {
        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData() << ": None\n";
    } 
    else {
        std::vector<std::string> path;
        for (int at = i_dest; at != -1; at = predecessors[at]) {
            path.push_back(vertices[at].getData());
        }

        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData() << ": ";
        for (int i = path.size() - 1; i >= 0; i--) {
            std::cout << path[i];
            if (i != 0) std::cout << " -> ";
        }
        std::cout << ". The length is " << distances[i_dest] << ". The cost is " << costs[i_dest] << ".\n";
    }
}

//Task 3: Find all shortest paths between the given origin airport and all the airports in the destination state. Dijkstra's Algorithm
template<typename T>
void Graph<T>::shortest_paths_to_state(const Vertex<T>& src, const std::string& dest_state, const std::vector<std::pair<std::string, std::string> >& airport_to_state) {
    int i_src = get_vertex_index(src);

    //If the origin airport does not exists, return
    if (i_src == -1) {
        std::cout << "Shortest route from " << src.getData() << " to " << dest_state << ": None\n";
        return;
    }

    //Marked all vertices as unvisited
    clean_visited();

    //Initialize
    std::vector<int> distances(vertices.size(), INT_MAX); //shortest distances
    std::vector<int> costs(vertices.size(), INT_MAX); //shortest costs
    std::vector<int> predecessors(vertices.size(), -1); 


    distances[i_src] = 0;
    costs[i_src] = 0;
    
    //Manually compare each neighbor visit distance, replace if found a better path
    MinHeap<Edge> heap;

    //start from the source
    vertices[i_src].setVisited(true);

    // Insert all initial edges from the source vertex
    for (const Edge& e : edges[i_src]) {
        heap.insert(e);
    }

    //While have not traverse through all edges
    while (!heap.is_empty()) {
        //Get the smallest cost
        Edge current = heap.delete_min();
        int u = current.src;
        int v = current.dest;

        //If destination already visited, skip
        if (vertices[v].getVisited()) continue;

        //Otherwise, update the infomation
        vertices[v].setVisited(true);
        distances[v] = distances[u] + current.distance;
        costs[v] = costs[u] + current.cost;
        predecessors[v] = u;

        // Add outgoing edges of vertex v to the heap
        for (const Edge& e : edges[v]) {
            if (!vertices[e.dest].getVisited()) {
                heap.insert(e);
            }
        }
    }

    // Output results for destination airports in the target state
    bool found = false;
    std::cout << "Origin -> Destination : Distance, Cost\n";

    for (int i = 0; i < vertices.size(); ++i) {
        if (distances[i] == INT_MAX) continue; // Skip unreachable vertices

        std::string airport_code = vertices[i].getData();

        // Match airport code with corresponding state
        for (const auto& pair : airport_to_state) {
            if (pair.first == airport_code && pair.second == dest_state) {
                found = true;

                // Reconstruct path from source to vertex
                std::vector<std::string> path;
                for (int at = i; at != -1; at = predecessors[at]) {
                    path.push_back(vertices[at].getData());
                }

                // Print path from source to dest
                for (int j = path.size() - 1; j >= 0; --j) {
                    std::cout << path[j];
                    if (j != 0) std::cout << " -> ";
                }

                // Output distance and cost of path
                std::cout << " : " << distances[i] << ", " << costs[i] << "\n";
                break;
            }
        }
    }

    if (!found) { // If no valid airport, print no paths
        std::cout << "No paths from " << src.getData() << " to state " << dest_state << " found.\n";
    }

    clean_visited();
    
}

// Task 4: Computes and displays the path between two airports
// within a certain number of stops.
template <typename T>
void Graph<T>::shortest_path_with_stops(const Vertex<T>& src, const Vertex<T>& dest, int n_stops) {
    int i_src = get_vertex_index(src); // Get index of source
    int i_dest = get_vertex_index(dest); // Get index of dest

    if (i_src == -1 || i_dest == -1) { // If either vertex is not in the graph, print no routes
        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData()
                  << " with max " << n_stops << " stops: None\n";
        return;
    }

    // Struct with all information needed for a path
    struct State {
        int node; // Current node index
        int dist; // Total distance
        int cost; // Total cost
        int stops; // Number of stops
        std::vector<int> path; // Path taken

        bool operator<(const State& other) const {
            return dist > other.dist; // Reversed for MinHeap
        }
    };

    MinHeap<State> heap;
    heap.insert({i_src, 0, 0, -1, {i_src}});

    int best_dist = INT_MAX;
    int best_cost = INT_MAX;
    std::vector<int> best_path; // Start from source with -1 stops

    while (!heap.is_empty()) { // Until the heap is empty,
        State current = heap.delete_min(); // Get state with shortest distance
    
        if (current.node == i_dest) { // If destination is found
            if (current.stops == n_stops) { // If number of stops is matched
                if (current.dist < best_dist) {
                    best_dist = current.dist;
                    best_cost = current.cost;
                    best_path = current.path;
                }
            }
            continue;
        }
    
        if (current.stops >= n_stops) { // Continue if path exceeds allowed stops
            continue;
        }
    
        for (const Edge& edge : edges[current.node]) { // Iterate edges for current node
            State next;
            next.node = edge.dest;
            next.dist = current.dist + edge.distance;
            next.cost = current.cost + edge.cost;
            next.stops = current.stops + 1;
            next.path = current.path;
            next.path.push_back(edge.dest);
    
            heap.insert(next);
        }
    }
    
    if (best_path.empty()) { // If no path was found, print none
        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData() << " with " << n_stops << " stops: None\n";
    } else { // Print the best path
        std::cout << "Shortest route from " << src.getData() << " to " << dest.getData() << " with " << n_stops << " stops: ";
        for (size_t i = 0; i < best_path.size(); i++) {
            std::cout << vertices[best_path[i]].getData();
            if (i != best_path.size() - 1) std::cout << " -> ";
        }
        // Print total length and cost of path
        std::cout << ". The length is " << best_dist << ". The cost is " << best_cost << ".\n";
    }
}

// Task 5: Counts the number of direct connections to each airport vertex,
// sorts them by connection, and prints them.
template <typename T>
void Graph<T>::count_direct_connections() {
    std::vector<std::pair<int, std::string>> connections; // Holds # of direct flights

    for (unsigned i = 0; i < vertices.size(); i++) { // Iterate through each airport vertex
        connections.push_back(std::make_pair(edges[i].size(), vertices[i].getData()));
    }

    std::vector<std::pair<int, std::string>> sorted_connections; // Holds # of connections sorted by connections
    int max = 0;
    auto max_i = connections.begin();
    while(!connections.empty()) {
        for (auto i = connections.begin(); i != connections.end(); i++) {
            if (i->first > max) { // Check if current connections is greater than max
                max = i->first; // If so, update max
                max_i = i; // Update iterator to max index
            }
        }

        sorted_connections.push_back(*max_i); // Push max element to sorted
        connections.erase(max_i);  // Remove max element from unsorted
        max = 0;
    }
    
    std::cout << "Airport  Connections" << std::endl;
    for (const auto& connection : sorted_connections) { // Iterate over each pair of connections and airport
        std::cout << "  " << connection.second << std::right << std::setw(11) << connection.first << std::endl; // Display airport code and connections
    }
}

// Task 6: Computes and displays an undirected graph
template<typename T>
Graph<T> Graph<T>::create_undirected_graph() {
    Graph<T> undirected; // Undirected graph

    for (const auto& vertex : this->vertices) { // Iterate through vertices in current graph
        undirected.insert_vertex(vertex); // Insert vertices into undirected graph
    }

    int v_size = vertices.size();

    // Track whether an edge was read
    std::vector<std::vector<bool>> processed(v_size, std::vector<bool>(v_size, false));

    for (int i = 0; i < v_size; i++) { // Iterate for each vertex i
        for (const auto& edge : edges[i]) { // Iterate through each vertex j
            int j = edge.dest; // Get destination vertex

            if (processed[i][j] || processed[j][i])
                continue; // Skip vertex if read already

            MinHeap<Edge> heap;
            heap.insert(edge); // Insert edge from i to j

            for (const auto& edge : edges[j]) { // Iterate through edges
                if (edge.dest == i) {
                    heap.insert(edge); // Insert opposite edge into heap
                    break;
                }
            }

            Edge minEdge = heap.delete_min(); // Get edge with the lowest cost

            // Add undirected edge both directions
            undirected.add_edge(vertices[i], vertices[j], minEdge.distance, minEdge.cost);
            undirected.add_edge(vertices[j], vertices[i], minEdge.distance, minEdge.cost);

            // Mark vertices as read
            processed[i][j] = true;
            processed[j][i] = true;
        }
    }

    return undirected;
}

// Task 7: Computes and displays the MST using Prim's algorithm
template<typename T>
void Graph<T>::prim_mst() {
    int n = vertices.size();
    //Check if G_u has vertices
    if (n == 0) {
        std::cout << "Graph is empty.\n";
        return;
    }

    // Set all vertices as unvisited
    clean_visited();

    //Initalization for MST
    std::vector<Edge> mstEdges; //vector to store all the edges with the smallest cost
    int totalCost = 0;
    MinHeap<Edge> edgeHeap;
    vertices[0].setVisited(true); //start from vertices[0]


    //Insert all edges of vertext [0] into the heap
    for (const Edge& e : edges[0]) {
        edgeHeap.insert(e);
    }

    // Perform Prim's algorithm

    //While we have not empty all edgeHeap and have not over flow mstEdge
    while (!edgeHeap.is_empty() && mstEdges.size() < n - 1) {
        //Find the egde with the smallest cost
        Edge minEdge;
        bool found = false;

        // Temporary vector store edges to reinsert later
        std::vector<Edge> temp;

        // Go through the edgeHeap
        while (!edgeHeap.is_empty()) {
            Edge e = edgeHeap.delete_min();

            //If the end of a vertice is not visited, check if the cost is smaller than the current cost in minEdge
            if (!vertices[e.dest].getVisited()) {
                if (!found || e.cost < minEdge.cost) {
                    if (found){
                        temp.push_back(minEdge);  // Save previous minimum edge
                    } 
                    //replace the previous value
                    minEdge = e;
                    found = true;
                } 
                else {
                    temp.push_back(e); //if not smaller, store in the temp vector for later use
                }
            }
        }

        // If no valid edge is found after traversed through all of edgeHeap (graph may be disconnected)
        if (!found) {
            std::cout << "No more edges to process or MST could not be constructed.\n";
            break;
        }

        // Reinsert remaining edges into the heap
        for (const Edge& e : temp) {
            edgeHeap.insert(e);
        }

        // Add the minEdge found earlier to mstEdge. Calculate the total cost
        vertices[minEdge.dest].setVisited(true);
        mstEdges.push_back(minEdge);
        totalCost += minEdge.cost;

        // Add all edges of the newly visited vertex to the heap
        for (const Edge& e : edges[minEdge.dest]) {
            if (!vertices[e.dest].getVisited()) {
                edgeHeap.insert(e);
            }
        }
    }

    //Sort MST edges by cost 
    for (int i = 0; i < mstEdges.size(); ++i) {
        int minIndex = i;
        for (int j = i + 1; j < mstEdges.size(); ++j) {
            if (mstEdges[j].cost < mstEdges[minIndex].cost) {
                minIndex = j;
            }
        }
        if (minIndex != i) {
            // Swap edges at i and minIndex
            Edge temp = mstEdges[i];
            mstEdges[i] = mstEdges[minIndex];
            mstEdges[minIndex] = temp;
        }
    }
    
    // If MST is not connected, exit
    if (mstEdges.empty()) {
        std::cout << "No MST found; graph may be disconnected or empty." << std::endl;
        return;
    }

    // Print the MST edges and total cost
    std::cout << "\nPrim's MST Edges:\n";
    for (const Edge& e : mstEdges) {
        std::cout << vertices[e.src].getData() << " -> " << vertices[e.dest].getData()
                  << " with cost " << e.cost << "\n";
    }
    std::cout << "Total cost of the MST: " << totalCost << "\n";
}

// Task 8: Computes and displays the MST using Kruskal's algorithm
template<typename T>
void Graph<T>::kruskal_mst() {
    clean_visited();

    std::vector<Edge> mstEdges; // Stores the edges for the MST
    std::vector<Edge> graphEdges; // Stores the edges from the graph
    int totalCost = 0;

    for (int u = 0; u < vertices.size(); u++) { // Iterate through vertices
        for (const Edge& edge : edges[u]) { // Iterate through edges
            int v = edge.dest;
            if (!vertices[v].getVisited()) {
                graphEdges.push_back(edge); // Add edge if vertex hasn't been visited
            }
        }
        vertices[u].setVisited(true);
    }

    DisjointSet ds(vertices.size());
    int edgesUsed = 0;
    
    while (!graphEdges.empty() && edgesUsed < vertices.size() - 1) { // While there are unchecked edges
        int minIndex = 0;
        // Find edge with minimum cost
        for (int i = 1; i < graphEdges.size(); ++i) { // Iterate through edges
            if (graphEdges[i].cost < graphEdges[minIndex].cost)
                minIndex = i; // Set index at element with lowest cost
        }

        // Remove edge with minimum cost
        Edge minEdge = graphEdges[minIndex];
        graphEdges.erase(graphEdges.begin() + minIndex);

        if (!ds.connected(minEdge.src, minEdge.dest)) { // If the edge connects a source and dest, add to MST
            ds.union_sets(minEdge.src, minEdge.dest); // Union source and dest
            mstEdges.push_back(minEdge); // Add edge to MST
            totalCost += minEdge.cost; // Add edge cost to total
            edgesUsed++; // Increment used edge count
        }
    }

    if (mstEdges.empty()) {
        std::cout << "No MST found; graph may be disconnected or empty." << std::endl;
        return;
    }

    // Print MST by printing all edges
    std::cout << "\nOrigin -> Destination: Cost:" << std::endl;
    for (const Edge& edge : mstEdges) {
        std::cout << vertices[edge.src].getData() << " -> " << vertices[edge.dest].getData() << ": " << edge.cost << std::endl;
    }
    std::cout << "\nTotal cost of MST: " << totalCost << std::endl;
}

