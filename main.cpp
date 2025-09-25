#include "Graph.h"
#include "Vertex.h"
#include "Graph.cpp"

#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;

// Finds and returns the index of an airport code vector if it exists
int find_index(const vector<Vertex<string>>& vertices, const string& airport_code) {
    for (unsigned i = 0; i < vertices.size(); i++) { // For each vertex
        if (vertices[i].getData() == airport_code) {
            return i; // Return index if airport has a vertex
        }
    }
    return -1; // Return -1 to indicate vertex does not exist
}

// Builds the directed weighted graph
tuple<Graph<string>, vector<Vertex<string>>, vector<pair<string, string>>> createDirectedGraph() {
    ifstream file("airports.txt");
    if (!file.is_open()) {
        cout << "Error opening file." << endl;
        return {};
    }

    Graph<string> G;
    vector<Vertex<string>> vertices;
    vector<pair<string, string>> airport_to_state;

    string origin_code, dest_code, origin_state, dest_state, num;
    int distance, cost, origin_index, dest_index;
    
    string line;
    getline(file, line);

    while (getline(file, line)) {
        stringstream ss(line);

        getline(ss, origin_code, ','); // Get origin airport code
        getline(ss, dest_code, ','); // Get dest airport code

        getline(ss, origin_state, ','); // Throw out city and comma
        getline(ss, origin_state, ' '); // Throw out space
        getline(ss, origin_state, '"'); // Get origin state

        getline(ss, dest_state, ','); // Throw out all up until comma
        getline(ss, dest_state, ','); // Throw out all up until comma
        getline(ss, dest_state, ' '); // Throw out space
        getline(ss, dest_state, '"'); // Get dest state

        getline(ss, num, ','); // Throw out comma
        getline(ss, num, ','); // Get cost
        distance = stoi(num);

        getline(ss, num, ','); // Throw out comma
        getline(ss, num, ','); // Get cost
        cost = stoi(num);

        origin_index = find_index(vertices, origin_code); // Find if vertex for airport exists
        if (origin_index == -1) { // If the origin vertex does not exist,
            vertices.emplace_back(origin_code); // Add vertex to vertices
            G.insert_vertex(vertices.back()); // Add vertex to graph
            origin_index = vertices.size() - 1; // Update position in vector
        }

        dest_index = find_index(vertices, dest_code); // Find if vertex for airport exists
        if (dest_index == -1) { // If the dest vertex does not exist,
            vertices.emplace_back(dest_code); // Add vertex to vertices
            G.insert_vertex(vertices.back()); // Add vertex to graph
            dest_index = vertices.size() - 1; // Update position in vector
        }

        G.add_edge(vertices[origin_index], vertices[dest_index], distance, cost); // Creates an edge between airport vertices

        bool origin_found = false, dest_found = false;
        
        for (auto &entry : airport_to_state) { // Iterate over each pair in the vector
            if (entry.first == origin_code) origin_found = true; // If the first entry is the origin code, mark as found
            if (entry.first == dest_code) dest_found = true; // If the first entry is the dest code, mark as found
        }

        // If the origin airport is not found, add to vector
        if (!origin_found) { airport_to_state.emplace_back(origin_code, origin_state); }
        // If the dest airport is not found, add to vector
        if (!dest_found) { airport_to_state.emplace_back(dest_code, dest_state); }
    }

    file.close();

    return make_tuple(G, vertices, airport_to_state);
}

// Task 1
void displayDirectedGraph(Graph<string> G) {
    cout << "Displaying weighted directed graph...\n" << endl;
    G.print();
}

// Task 2
void displayShortestPath(Graph<string>& G) {
    string origin, dest;
    cout << "Enter the origin airport code: ";
    cin >> origin;
    cout << "Enter the destination airport code: ";
    cin >> dest;
    cout << "Finding shortest path from " << origin << " to " << dest << "...\n" << endl;
    
    G.shortest_path(origin, dest);
}

// Task 3
void displayShortestPathsToState(Graph<string>& G, vector<pair<string, string>>& airport_to_state) {
    string origin, state;
    cout << "Enter the origin airport code: ";
    cin >> origin;
    cout << "Enter the destination state (abbreviation): ";
    cin >> state;
    cout << "Finding shortest paths from " << origin << " to airports in " << state << "...\n" << endl;
    
    G.shortest_paths_to_state(origin, state, airport_to_state);
}

// Task 4
void displayShortestPathWithStops(Graph<string>& G) {
    string origin, dest;
    int stops;
    cout << "Enter the origin airport code: ";
    cin >> origin;
    cout << "Enter the destination airport code: ";
    cin >> dest;
    cout << "Enter the maximum number of stops: ";
    cin >> stops;
    cout << "Finding shortest paths from " << origin << " to " << dest << " in " << stops << " stops...\n" << endl;
    
    G.shortest_path_with_stops(origin, dest, stops);
}

// Task 5
void displayDirectConnections(Graph<string> G) {
    cout << "Displaying total direct flights connections to each airport..." << endl;
    G.count_direct_connections();
}

// Task 6
void displayUndirectedGraph(Graph<string>& G_u) {
    cout << "Displaying undirected graph..." << endl;
    G_u.print();
}

// Task 7
void displayPrimMST(Graph<string>& G_u) {
    cout << "Displaying MST generated by Prim's algorithm..." << endl;
    G_u.prim_mst();
}

// Task 8
void displayKruskalMST(Graph<string>& G_u) {
    cout << "Displaying MST generated by Kruskal's algorithm..." << endl;
    G_u.kruskal_mst();
}

int main() {
    int option;
    bool flag = true;

    auto result = createDirectedGraph();
    Graph<string> G = get<0>(result);
    vector<Vertex<string>> vertices = get<1>(result);
    vector<pair<string, string>> airport_to_state = get<2>(result);
    
    Graph<string> G_u = G.create_undirected_graph();

    while(flag) {
        cout << "\n--- Flight Route Optimizer ---" << endl;
        cout << "[1] Display weighted directed graph" << endl;
        cout << "[2] Display the shortest path between two airports" << endl;
        cout << "[3] Display all shortest paths between an airport and all airports in a state" << endl;
        cout << "[4] Display the shortest path between two airports in a given number of stops" << endl;
        cout << "[5] Display the total direct flight connections to each airport" << endl;
        cout << "[6] Display undirected graph" << endl;
        cout << "[7] Display an MST generated by Prim's algorithm" << endl;
        cout << "[8] Display an MST generated by Kruskal's algorithm" << endl;
        cout << "[0] Quit" << endl;
        cout << "Enter an option: ";
        cin >> option;
        cout << endl;

        switch(option) {
            case 1: displayDirectedGraph(G); break;
            case 2: displayShortestPath(G); break;
            case 3: displayShortestPathsToState(G, airport_to_state); break;
            case 4: displayShortestPathWithStops(G); break;
            case 5: displayDirectConnections(G); break;
            case 6: displayUndirectedGraph(G_u); break;
            case 7: displayPrimMST(G_u); break;
            case 8: displayKruskalMST(G_u); break;
            case 0: cout << "Exiting program..." << endl; flag = false; break;
            default: cout << "Invalid option. Please try again" << endl;
        }
    }

    return 0;
}
