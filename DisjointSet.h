class DisjointSet {
    std::vector<int> parent; // Stores parent of each element

public:
    DisjointSet(int size) {
        parent.resize(size);
        for (int i = 0; i < size; ++i)
            parent[i] = i; // Set each element as its own parent
    }

    // Finds the root of x and compresses path
    int find(int x) {
        if (parent[x] != x) // If x is not the root
            parent[x] = find(parent[x]); // Point x to root
        return parent[x]; // Return root
    }

    // Merges sets containing x and y
    void union_sets(int x, int y) {
        int root_x = find(x); // Find root of x
        int root_y = find(y); // Find root of y
        if (root_x != root_y) // If x and y are in different sets
            parent[root_y] = root_x; // Merge y roo to x root
    }

    bool connected(int x, int y) {
        return find(x) == find(y); // Return true if both elements have the same root
    }
};
