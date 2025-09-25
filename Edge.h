#ifndef EDGE_H
#define EDGE_H

class Edge {
public:
    Edge(int src = 0, int dest = 0, int distance = 0, int cost = 0) : 
        src(src), dest(dest), distance(distance), cost(cost) {}

    bool operator<(const Edge& other) const {
        if(this->distance == other.distance) {
            return cost < other.cost;
        }

        return distance < other.distance;
    }

    int src;
    int dest;
    int distance;
    int cost;
};

#endif