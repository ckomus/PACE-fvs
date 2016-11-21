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

#ifndef GRAPH_HH_INCLUDED
#define GRAPH_HH_INCLUDED

#include <cassert>

#include <set>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <vector>

typedef unsigned Vertex;
constexpr Vertex NULL_VERTEX = -1;

struct Edge {
    Edge() { }
    Edge(Vertex nu, Vertex nv) {
	if (nu < nv) {
	    u = nu;
	    v = nv;
	} else {
	    u = nv;
	    v = nu;
	}
    }
    bool operator==(const Edge& other) const { return u == other.u && v == other.v; }
    bool operator<(const Edge& other) const { if (u != other.u) return u < other.u; else return v < other.v; }

    Vertex u, v;
};


class Graph {
public:
    typedef std::multiset<Vertex> Neighbors;
    typedef std::map<Vertex, Neighbors> NeighborMap;

    // constructors
    Graph() {}
    Graph(std::string g6);
    static std::pair<Graph, std::map<Vertex, std::string> > read(std::istream& in);

    // query
    std::size_t numVertices() const { return neighbors_.size(); }
    size_t numEdges() const {
        size_t count = 0;
        for (auto it : neighbors_) {
            count += it.second.size();
	    count += it.second.count(it.first); // also double-count self-loops
	}
        assert(count % 2 == 0);
        return count / 2;
    }
    bool hasVertex(Vertex u) const { return neighbors_.find(u) != neighbors_.end(); }
    size_t deg(Vertex u) const { return neighbors(u).size(); }
    bool hasEdge(Vertex u, Vertex v) const {
        NeighborMap::const_iterator nu = neighbors_.find(u);
        if (nu == neighbors_.end())
            throw std::out_of_range("Graph::hasEdge " + std::to_string(u));
        NeighborMap::const_iterator nv = neighbors_.find(v);
        if (nv == neighbors_.end())
            throw std::out_of_range("Graph::hasEdge " + std::to_string(v));
        return nu->second.find(v) != nu->second.end();
    }
    bool hasEdge(Edge e) const { return hasEdge(e.u, e.v); }
    const NeighborMap& neighbors() const { return neighbors_; }
    const Neighbors& neighbors(Vertex u) const {
        NeighborMap::const_iterator it = neighbors_.find(u);
        if (it == neighbors_.end())
            throw std::out_of_range("Graph::neighbors " + std::to_string(u));
        return it->second;
    }

    // iterators
    class Vertices {
    public:
        Vertices(const Graph& g) : g_(g) { }
        class Iterator {
        public:
            Iterator(NeighborMap::const_iterator it) : it_(it) { }
            bool operator!=(const Iterator& other) const { return it_ != other.it_; }
            Vertex operator*() const { return it_->first; }
            const Iterator& operator++() { ++it_; return *this; }
        private:
            NeighborMap::const_iterator it_;
        };
        Iterator begin() const { return Iterator(g_.neighbors().begin()); }
        Iterator end() const { return Iterator(g_.neighbors().end()); }
    private:
        const Graph& g_;
    };
    Vertices vertices() const { return Vertices(*this); }

    class VertexPairs {
    public:
	VertexPairs(const Graph& g) : g_(g) { }
	class Iterator {
	public:
	    Iterator(NeighborMap::const_iterator begin, NeighborMap::const_iterator end)
		: it1_(begin), it2_(begin), end_(end) {
		if (it1_ != end_)
		    if (++it2_ == end)
			++it1_;
	    }
	    // for comparing against end() only
	    bool operator!=(const Iterator& other) const { return it1_ != other.it1_; }
	    Edge operator*() const { return {it1_->first, it2_->first}; }
	    const Iterator& operator++() {
		while (++it2_ == end_) {
		    ++it1_;
		    it2_ = it1_;
		    if (it2_ == end_)
			break;
		}
		return *this;
	    }
	private:
	    NeighborMap::const_iterator it1_, it2_, end_;
	};
	Iterator begin() const { return Iterator(g_.neighbors().begin(), g_.neighbors().end()); }
	Iterator end() const { return Iterator(g_.neighbors().end(), g_.neighbors().end()); }
    private:
	const Graph& g_;
    };
    VertexPairs vertexPairs() const { return VertexPairs(*this); }

    class Edges {
    public:
        Edges(const Graph& g) : g_(g) { }
        class Iterator {
        public:
            Iterator(NeighborMap::const_iterator it, NeighborMap::const_iterator it_end)
                : it_(it), it_end_(it_end) {
                while (it_ != it_end_) {
                    nit_ = it_->second.lower_bound(it_->first);
                    if (nit_ != it_->second.end())
                        break;
                    ++it_;
                }
            }
            bool operator!=(const Iterator& /*other*/) const {
                // only for comparison against end()
                return it_ != it_end_;
            }
            Edge operator*() const { return {it_->first, *nit_}; }
            const Iterator& operator++() {
                ++nit_;
                while (nit_ == it_->second.end()) {
                    ++it_;
                    if (it_ == it_end_)
                        break;
                    nit_ = it_->second.lower_bound(it_->first);
                }
                return *this;
            }
        private:
            NeighborMap::const_iterator it_, it_end_;
            Neighbors::const_iterator nit_;
        };

        Iterator begin() const { return Iterator(g_.neighbors().begin(), g_.neighbors().end()); }
        Iterator end() const { return Iterator(g_.neighbors().end(),  g_.neighbors().end()); }

    private:
        const Graph& g_;
    };
    Edges edges() const { return Edges(*this); }

    // modification
    Vertex addVertex() {
        Vertex u;
        if (neighbors_.empty())
            u = 0;
        else
            u = neighbors_.rbegin()->first + 1;
        neighbors_.insert(neighbors_.end(), {u, Neighbors()});
        return u;
    }

    void addEdge(Vertex u, Vertex v) {
        auto nu_it = neighbors_.find(u);
        if (nu_it == neighbors_.end())
            throw std::out_of_range("Graph::addEdge " + std::to_string(u));
        Neighbors::iterator lb = nu_it->second.lower_bound(v);
        nu_it->second.insert(lb, v);

	if (u != v) {
	    auto nv_it = neighbors_.find(v);
	    if (nv_it == neighbors_.end())
		throw std::out_of_range("Graph::addEdge " + std::to_string(v));
	    nv_it->second.insert(u);
	}
    }

    void deleteEdge(Vertex u, Vertex v) {
        auto nu_it = neighbors_.find(u);
        if (nu_it == neighbors_.end())
            throw std::out_of_range("Graph::deleteEdge " + std::to_string(u));
        Neighbors::iterator pv = nu_it->second.find(v);
        if (pv == nu_it->second.end())
            throw std::invalid_argument("Graph::deleteEdge: not an edge {"
                                        + std::to_string(u) + ", " +  std::to_string(v) + '}');
        nu_it->second.erase(pv);

	if (u != v) {
	    auto nv_it = neighbors_.find(v);
	    assert(nv_it != neighbors_.end());
	    Neighbors::iterator pu = nv_it->second.find(u);
	    assert(pu != nv_it->second.end());
	    nv_it->second.erase(pu);
	}
    }

    void addEdge   (const Edge& e) { addEdge   (e.u, e.v); }
    void deleteEdge(const Edge& e) { deleteEdge(e.u, e.v); }

    void deleteVertex(Vertex u) {
	for (Vertex v : neighbors(u))
	    if (v != u)
		neighbors_[v].erase(u);
	neighbors_.erase(u);
    }

    bool isConnected() const;
    std::vector<Graph> connectedComponents() const;
    std::vector<Edge> bridges() const;
    Graph subgraph(const std::set<Vertex> vs) const;

private:
    NeighborMap neighbors_;
};

#endif // GRAPH_HH_INCLUDED
