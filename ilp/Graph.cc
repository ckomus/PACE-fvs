/* commode -- find significant protein complex graph models
   Copyright (C) 2013, 2014, 2015  Falk Hüffner

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.  */

#include <algorithm>

#include "Graph.hh"
#include "Util.hh"

Graph::Graph(std::string g6) {
    for (size_t i = 0; i < g6.size(); ++i)
	g6[i] -= 63;
    size_t n = g6[0];
    for (size_t i = 0; i < n; ++i) {
	Vertex u = addVertex();
	assert(i == u);
    }

    size_t b = 0;
    for (size_t j = 0 ; j < n; ++j) {
	for (size_t i = 0; i < j; ++i) {
	    size_t byte = 1 + (b / 6);
	    size_t bit = 5 - (b % 6);
	    if ((g6[byte] >> bit) & 1)
		addEdge(i, j);
	    ++b;
	}
    }
};

std::pair<Graph, std::map<Vertex, std::string>> Graph::read(std::istream& in) {
    Graph g;
    std::map<std::string, Vertex> vertexNumbers;
    std::map<Vertex, std::string> vertexNames;
    std::string line;
    std::size_t linenr = 0;
    while (std::getline(in, line)) {
	++linenr;
	std::string::size_type hash = line.find_first_of("#%");
	if (hash != std::string::npos)
	    line.resize(hash);
	std::vector<std::string> fields = Util::split(line);
	if (fields.size() == 0)
	    continue;
	if (fields.size() != 2)
	    throw std::runtime_error("line " + std::to_string(linenr) + ": "
				     + std::to_string(fields.size()) + " fields, but 2 required");
	const std::string& nu = fields[0], nv = fields[1];
	Vertex u, v;
	auto pu = vertexNumbers.find(nu);
	if (pu != vertexNumbers.end())
	    u = pu->second;
	else {
	    vertexNumbers[nu] = u = g.addVertex();
	    vertexNames[u] = nu;
	}
	auto pv = vertexNumbers.find(nv);
	if (pv != vertexNumbers.end())
	    v = pv->second;
	else {
	    vertexNumbers[nv] = v = g.addVertex();
	    vertexNames[v] = nv;
	}
	g.addEdge(u, v);
    }
    return {g, vertexNames};
}

bool Graph::isConnected() const {
    if (numVertices() == 0)
        return true;
    std::set<Vertex> seen;
    Vertex start = neighbors_.begin()->first;
    seen.insert(start);
    std::vector<Vertex> queue;
    queue.push_back(start);
    while (!queue.empty()) {
        Vertex u = queue.back();
        queue.pop_back();
        for (Vertex v : neighbors(u)) {
            if (seen.find(v) == seen.end()) {
                seen.insert(v);
                queue.push_back(v);
            }
        }
    }
    return seen.size() == numVertices();
}

void dfs(const Graph& g,
	 std::map<Vertex, int>& dfsNumber,
	 std::map<Vertex, Vertex>& dfsParent,
	 std::map<Vertex, int>& dfsOrder,
	 std::map<Vertex, std::set<Vertex>>& backEdges,
	 int u, int parent, int& d) {
    dfsOrder[d] = u;
    dfsNumber[u] = ++d;
    dfsParent[u] = parent;
    for (int v : g.neighbors(u)) {
	if (v == parent)
	    continue;
	if (dfsNumber.find(v) != dfsNumber.end())
	    backEdges[v].insert(u);
	else
	    dfs(g, dfsNumber, dfsParent, dfsOrder, backEdges, v, u, d);
    }
}

std::vector<Graph> Graph::connectedComponents() const {
    std::set<Vertex> seen;
    std::vector<Graph> ccs;
    for (Vertex u : vertices()) {
        if (seen.find(u) != seen.end())
            continue;
        seen.insert(u);
        std::set<Vertex> cc;
        std::vector<Vertex> queue;
        queue.push_back(u);
        cc.insert(u);
        while (!queue.empty()) {
            Vertex u = queue.back();
            queue.pop_back();
            for (Vertex v : neighbors(u)) {
                if (seen.find(v) == seen.end()) {
                    seen.insert(v);
                    queue.push_back(v);
                    cc.insert(v);
                }
            }
        }
	ccs.push_back(subgraph(cc));
    }
    return ccs;
}

std::vector<Edge> Graph::bridges() const {
    int n = numVertices();
    if (n < 2)
	return {};
    std::map<Vertex, int> dfsNumber;
    std::map<Vertex, Vertex> dfsParent;
    std::map<Vertex, std::set<Vertex>> backEdges;
    std::map<Vertex, int> dfsOrder;
    int d = 0;
    for (Vertex u : vertices()) {
	if (dfsNumber.find(u) == dfsNumber.end())
	    dfs(*this, dfsNumber, dfsParent, dfsOrder, backEdges, u, u, d);
    }
    assert(d == n);
    std::set<Vertex> visited;
    std::set<Edge> traversed;
    for (int i = 0; i < n; ++i) {
	Vertex u = dfsOrder[i];
	if (!backEdges[u].empty()) {
	    visited.insert(u);
	    for (Vertex v : backEdges[u]) {
		traversed.insert({u, v});
		Vertex w = v;
		while (visited.find(w) == visited.end()) {
		    visited.insert(w);
		    traversed.insert({w, dfsParent[w]});
		    w = dfsParent[w];
		}
	    }
	}
    }
    std::vector<Edge> bridges;
    for (Edge e : edges())
	if (traversed.find(e) == traversed.end())
	    bridges.push_back(e);
    return bridges;
}

template<class InputIt1, class InputIt2, class OutputIt>
OutputIt set_multiset_intersection(InputIt1 set, InputIt1 set_end,
				   InputIt2 multiset, InputIt2 multiset_end,
				   OutputIt out) {
    while (set != set_end && multiset != multiset_end) {
        if (*set < *multiset) {
            ++set;
        } else  {
            if (!(*multiset < *set))
                *out++ = *set;
            ++multiset;
        }
    }
    return out;
}

Graph Graph::subgraph(const std::set<Vertex> vs) const {
    Graph g;
    for (Vertex u : vs) {
	NeighborMap::iterator it = g.neighbors_.insert({u, {}}).first;
	const Neighbors& n_u = neighbors(u);
	set_multiset_intersection(vs.begin(), vs.end(),
				  n_u.begin(), n_u.end(),
				  std::inserter(it->second, it->second.end()));
    }
    return g;
}
