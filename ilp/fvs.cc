/* fvs -- find a Feedback Vertex Set
   Copyright (C) 2015, 2016  Falk Hüffner

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

#include "Graph.hh"
#include "LP/LP.hh"
#include "Util.hh"

#include <ctime>

#include <unistd.h>

bool verbose = false;
double timeLimit = 0;

void usage(const std::string& progname) {
    std::cerr << "usage: " << progname << " [-s] [-t S] [-v] < graph\n";
}

bool sim(double x, double y) {
    return fabs(x - y) < 1e-6;
}

double now() {
    return double(std::clock()) / CLOCKS_PER_SEC;
}

// O(n^2) additive +1 approximation for shortest cycle: for each
// vertex in the graph, start a BFS until the first cycle is closed
// (then stop and move on to the next vertex); return the shortest
// cycle found.
std::vector<Vertex> bfs(const Graph& g, Vertex s) {
    std::vector<Vertex> queue;
    queue.push_back(s);
    std::map<Vertex, Vertex> pred;
    pred[s] = NULL_VERTEX;
    while (!queue.empty()) {
	Vertex u = queue.back();
	queue.pop_back();
	for (Vertex v : g.neighbors(u)) {
	    if (pred.find(v) == pred.end()) {
		queue.push_back(v);
		pred[v] = u;
	    } else if (v != pred[u]) {
		std::vector<Vertex> cu, cv;
		for (Vertex w = u; w != NULL_VERTEX; w = pred[w])
		    cu.push_back(w);
		for (Vertex w = v; w != NULL_VERTEX; w = pred[w])
		    cv.push_back(w);
		while (cu.size() > 1 && cv.size() > 1 && cu[cu.size() - 2] == cv[cv.size() - 2]) {
		    cu.pop_back();
		    cv.pop_back();
		}
		cv.pop_back();
		for (Vertex u : cv)
		    cu.push_back(u);
		return cu;
	    }
	}
    }
    return {};
}

std::vector<std::vector<Vertex>> findCycles(const Graph& g) {
    std::vector<std::vector<Vertex>> shortestCycles;
    for (Vertex u : g.vertices()) {
	auto cycle = bfs(g, u);
	if (cycle.empty())
	    continue;
	if (shortestCycles.empty() || cycle.size() < shortestCycles[0].size())
	    shortestCycles = {cycle};
	else if (cycle.size() == shortestCycles[0].size())
	    shortestCycles.push_back(cycle);
    }
    return shortestCycles;
}

//////////////////////////////////////////////////////////////////////

typedef std::vector<Vertex> Vertices;
typedef std::vector<Vertices> VertexSets;

std::map<Vertex, int> degreeLabeling(const Graph& g) {
    constexpr bool white = false, black = true;
    std::map<Vertex, size_t> degree;
    std::map<Vertex, bool> color;
    std::map<Vertex, int> l;
    auto n = g.numVertices();
    for (auto v : g.vertices()) {
	color[v] = white;
	degree[v] = g.deg(v);
    }
    for (size_t i = 1; i <= n; ++i) {
	auto min_degree = n;
	Vertex v;
	for (auto x : g.vertices()) {
	    if (color[x] == white && degree[x] < min_degree) {
		v = x;
		min_degree = degree[x];
	    }
	}
	l[v] = i;
	color[v] = black;
	for (auto u : g.neighbors(v))
	    if (color[u] == white)
		--degree[u];
    }
    return l;
}

std::pair<VertexSets, VertexSets> triplets(const Graph& g, const std::map<Vertex, int>& l) {
    VertexSets t, c;
    for (Vertex u : g.vertices()) {
	for (Vertex x : g.neighbors(u)) {
	    for (Vertex y : g.neighbors(u)) {
		if (l.at(u) < l.at(x) && l.at(x) < l.at(y)) {
		    if (!g.hasEdge(x, y))
			t.push_back({x, u, y});
		    else
			c.push_back({x, u, y});
		}
	    }
	}
    }
    return {t, c};
}

void blockNeighbors(const Graph& g, std::map<Vertex, int>& blocked, Vertex v) {
    for (Vertex u : g.neighbors(v))
	++blocked[u];
}

void unblockNeighbors(const Graph& g, std::map<Vertex, int>& blocked, Vertex v) {
    for (Vertex u : g.neighbors(v))
	if (blocked[u] > 0)
	    --blocked[u];
}

void ccVisit(const Graph& g, const std::map<Vertex, int>& l,
	     const Vertices& p, VertexSets& c, int key, std::map<Vertex, int>& blocked) {
    auto u_t = p.back();
    blockNeighbors(g, blocked, u_t);
    for (auto v : g.neighbors(u_t)) {
	if (l.at(v) > key && blocked[v] == 1) {
	    Vertices p1 = p;
	    p1.push_back(v);
	    if (g.hasEdge(v, p[0]))
		c.push_back(p1);
	    else
		ccVisit(g, l, p1, c, key, blocked);
	}
    }
    unblockNeighbors(g, blocked, u_t);
}

VertexSets chordlessCycles(const Graph& g) {
    auto l = degreeLabeling(g);
    VertexSets t, c;
    std::tie(t, c) = triplets(g, l);
    std::map<Vertex, int> blocked;
    for (auto u : g.vertices())
	blocked[u] = 0;
    while (!t.empty()) {
	auto p = t.back();
	t.pop_back();
	auto u = p[1];
	blockNeighbors(g, blocked, u);
	ccVisit(g, l, p, c, l[u], blocked);
    }
    return c;
}

//////////////////////////////////////////////////////////////////////

std::pair<Graph, std::vector<Vertex>>
reduce(Graph g, const std::map<Vertex, std::string>& /*vertexNames*/) {
    std::vector<Vertex> fvs;
    //return {g, fvs};
    bool debugDR = true;
    while (true) {
	std::set<Vertex> todo(g.vertices().begin(), g.vertices().end());
	while (!todo.empty()) {
	    Vertex u = Util::SetOps::pop(todo);
	    if (!g.hasVertex(u))
		continue;
	    if (g.hasEdge(u, u)) {
		// Rule: if there is a self-loop {u, u}, then put u into the fvs
		fvs.push_back(u);
		for (Vertex v : g.neighbors(u))
		    todo.insert(v);
		g.deleteVertex(u);
		if (debugDR) std::cerr << "remove self-loop " << u << ' ' << u << std::endl;
		continue;
	    }

	    // Rule: if u has degree 0 or 1, then delete u
	    if (g.deg(u) <= 1) {
		for (Vertex v : g.neighbors(u))
		    todo.insert(v);
		if (debugDR) std::cerr << "remove degree-" << g.deg(u) << " vertex " << u << std::endl;
		g.deleteVertex(u);
		continue;
	    }
	    // Rule: shortcut degree-2 vertex
	    if (g.deg(u) == 2) {
		Vertex v = *g.neighbors(u).begin();
		Vertex w = *g.neighbors(u).rbegin();
		g.deleteVertex(u);
		g.addEdge(v, w);
		todo.insert(v);
		todo.insert(w);
		if (debugDR) std::cerr << "shortcut degree-2 vertex " << u << std::endl;
		continue;
	    }
	}

	// reduce multiplicity of edges to at most two
	std::map<Edge, size_t> edges;
	for (Edge e : g.edges())
	    ++edges[e];
	bool reduced = false;
	for (auto p : edges) {
	    while (p.second > 2) {
		g.deleteEdge(p.first);
		--p.second;
		reduced = true;
		if (debugDR) std::cerr << "reduce multiplicity " << p.first.u << ' ' << p.first.v << std::endl;
	    }
	}
	if (reduced)
	    continue;

	// delete bridges
	auto bridges = g.bridges();
	if (!bridges.empty()) {
	    for (Edge e : bridges) {
		g.deleteEdge(e);
		if (debugDR) std::cerr << "delete bridge " << e.u << ' ' << e.v << std::endl;
	    }
	    continue;
	}
	break;
    }
    return {g, fvs};
}

std::vector<Vertex> fvsHeuristic(const Graph& g0) {
    Graph g = g0;
    std::vector<Vertex> fvs;

    double t = now();
    while (true) {
       auto cycles = findCycles(g);
       if (cycles.empty())
	   return fvs;
       for (auto& cycle : cycles) {
	   size_t maxDeg = 0;
	   Vertex best = NULL_VERTEX;
	   bool doneAlready = false;
	   for (Vertex u : cycle) {
	       if (!g.hasVertex(u)) {
		   doneAlready = true;
		   break;
	       }
	       if (g.deg(u) > maxDeg) {
		   maxDeg = g.deg(u);
		   best = u;
	       }
	   }
	   if (doneAlready)
	       continue;
	   fvs.push_back(best);
	   g.deleteVertex(best);
	   if (verbose && now() - t > 2) {
	       std::cerr << "heuristic running, vertices left: " << g.numVertices()
			 << " fvs size: " << fvs.size() << std::endl;
	       t = now();
	   }
       }
    }
    return fvs;
}

size_t fvsLB(const Graph& g0) {
    Graph g = g0;
    size_t lb = 0;
    while (true) {
       auto cycles = findCycles(g);
       if (cycles.empty())
	   break;
       for (auto& cycle : cycles) {
	   bool doneAlready = false;
	   for (Vertex u : cycle) {
	       if (!g.hasVertex(u)) {
		   doneAlready = true;
		   break;
	       }
	   }
	   if (doneAlready)
	       continue;
	   for (Vertex u : cycle)
	       g.deleteVertex(u);
	   ++lb;
       }
    }
    return lb;
}

class Callback : public LP::Callback {
public:
    Callback(const Graph& g, const std::map<Vertex, LP::Var>& vars,
					      size_t nother_lb, size_t nother_ub,
					      size_t nbest_lb, size_t nbest_ub,
					      size_t nbest_comp_lb, size_t nbest_comp_ub)
        : g_(g), vars_(vars),
	  other_lb(nother_lb), other_ub(nother_ub),
	  best_lb(nbest_lb), best_ub(nbest_ub),
	  best_comp_lb(nbest_comp_lb), best_comp_ub(nbest_comp_ub) { }
    void call(LP::CallbackModel& model) {
#if 1
	double flb = model.bestBound();
	if (flb > 0) {
	    if (sim(flb, nearbyint(flb)))
		flb = nearbyint(flb);
	    else
		flb = ceil(flb);
	    size_t lb = flb;
	    if (lb > best_comp_lb) {
		std::cerr << "reduced-k-lb: " << lb << std::endl;
		best_comp_lb = lb;
	    }
	    if (other_lb + lb > best_lb) {
		std::cerr << "lb: " << other_lb + lb << std::endl;
		best_lb = other_lb + lb;
	    }
	}

#endif
        Graph g2 = g_;
        for (Vertex u : g_.vertices())
            if (model.value(vars_.at(u)))
                g2.deleteVertex(u);

	auto cycles = findCycles(g2);

	if (cycles.empty()) {
	    double fub = model.bestObjectiveValue();
	    if (fub < LP::inf) {
		assert(sim(fub, nearbyint(fub)));
		size_t ub = nearbyint(fub);
		if (ub < best_comp_ub) {
		    std::cerr << "reduced-k-ub: " << ub << std::endl;
		    best_comp_ub = ub;
		}
		if (other_ub + ub < best_ub) {
		    std::cerr << "ub: " << other_ub + ub << std::endl;
		    best_ub = other_ub + ub;
		}
	    }
	    return;
	}

	std::set<std::set<Vertex>> used;
	for (auto c : cycles) {
	    std::set<Vertex> vs(c.begin(), c.end());
	    if (used.find(vs) != used.end())
		continue;
	    used.insert(vs);
	    LP::Expr sum;
	    for (Vertex u : c)
		sum += vars_.at(u);
	    model.addLazyConstraint(sum >= 1);
	}
    }

private:
    const Graph& g_;
    const std::map<Vertex, LP::Var>& vars_;
    size_t other_lb, other_ub;
    size_t best_lb, best_ub;
    size_t best_comp_lb, best_comp_ub;
};

std::pair<size_t, std::vector<Vertex>> fvsHS(const Graph& g, const std::vector<Vertex>& heurFvs,
					     size_t /*other_lb*/, size_t /*other_ub*/,
					     size_t /*best_lb*/, size_t /*best_ub*/,
					     size_t /*best_comp_lb*/, size_t /*best_comp_ub*/,
					     LP::Solver solver = LP::Solver::GUROBI) {
    auto cycles = chordlessCycles(g);
    LP::Model& model = *LP::newModel(solver);
    std::map<Vertex, LP::Var> vars;
    LP::Expr objective;
    for (Vertex u : g.vertices()) {
	vars[u] = model.addBoolVar("v_" + std::to_string(u));
	objective += vars[u];
    }
    model.minimize(objective);

    for (auto cycle : cycles) {
	LP::Expr lhs;
	for (auto u : cycle)
	    lhs += vars[u];
	model.addConstraint(lhs >= 1);
    }

    if (verbose)
	model.setVerbose();

    if (!heurFvs.empty()) {
	std::map<LP::Var, double> mipStart;
	for (Vertex u : g.vertices())
	    mipStart[vars[u]] = 0;
	for (Vertex u : heurFvs)
	    mipStart[vars[u]] = 1;
	model.MIPstart(mipStart);
    }

    if (timeLimit)
	model.setTimeLimit(timeLimit);
    model.solve();
    size_t lb = 0;
    std::vector<Vertex> fvs;

    if (model.status() == LP::Status::OPTIMAL || model.status() == LP::Status::ABORTED_FEASIBLE) {
	for (Vertex u : g.vertices()) {
	    //std::cerr << "vars[" << u << "].value() = " << vars[u].value() << std::endl;
	    if (vars[u].value()) {
		fvs.push_back(u);
	    }
	}
	assert(sim(fvs.size(), model.objectiveValue()));
    }

    switch (model.status()) {
    case LP::Status::OPTIMAL:
	lb = fvs.size();
	break;
    case LP::Status::ABORTED_FEASIBLE: case LP::Status::ABORTED_NO_FEASIBLE: {
	double flb = model.bestBound();
	if (sim(flb, round(flb)))
	    lb = round(flb);
	else
	    lb = ceil(flb);
	lb = model.bestBound();
	if (model.status() == LP::Status::ABORTED_NO_FEASIBLE)
	    fvs = heurFvs;
	//std::cerr << "bound = " << model.bestBound()
	//	  << " lb = " << lb << " fvs.size() = " << fvs.size() << std::endl;
	assert(lb < fvs.size());
	break;
        }
    default:
	abort();
    }

    delete &model;
    return {0, fvs};
}

std::pair<size_t, std::vector<Vertex>> fvsILP(const Graph& g, const std::vector<Vertex>& heurFvs,
					      size_t other_lb, size_t other_ub,
					      size_t best_lb, size_t best_ub,
					      size_t best_comp_lb, size_t best_comp_ub,
					      LP::Solver solver = LP::Solver::GUROBI) {
    LP::Model& model = *LP::newModel(solver);
    // we are more interested in upper bounds than lower bounds
    model.setEmphasis(LP::Emphasis::FIND_FEASIBLE);
    std::map<Vertex, LP::Var> vars;
    LP::Expr objective;
    for (Vertex u : g.vertices()) {
	assert(g.deg(u) >= 3);
	vars[u] = model.addBoolVar("v_" + std::to_string(u));
	objective += vars[u];
    }
    model.minimize(objective);
    if (verbose)
	model.setVerbose();

    if (!heurFvs.empty()) {
	std::map<LP::Var, double> mipStart;
	for (Vertex u : g.vertices())
	    mipStart[vars[u]] = 0;
	for (Vertex u : heurFvs)
	    mipStart[vars[u]] = 1;
	model.MIPstart(mipStart);
    }

    Callback callback(g, vars, other_lb, other_ub, best_lb, best_ub, best_comp_lb, best_comp_ub);
    model.setCallback(callback, LP::LAZY_CONSTRAINTS);

    if (timeLimit)
	model.setTimeLimit(timeLimit);
    model.solve();
    size_t lb = 0;
    std::vector<Vertex> fvs;

    if (model.status() == LP::Status::OPTIMAL || model.status() == LP::Status::ABORTED_FEASIBLE) {
	for (Vertex u : g.vertices()) {
	    //std::cerr << "vars[" << u << "].value() = " << vars[u].value() << std::endl;
	    if (vars[u].value()) {
		fvs.push_back(u);
	    }
	}
	assert(sim(fvs.size(), model.objectiveValue()));
    }

    switch (model.status()) {
    case LP::Status::OPTIMAL:
	lb = fvs.size();
	break;
    case LP::Status::ABORTED_FEASIBLE: case LP::Status::ABORTED_NO_FEASIBLE: {
	double flb = model.bestBound();
	if (sim(flb, round(flb)))
	    lb = round(flb);
	else
	    lb = ceil(flb);
	lb = model.bestBound();
	if (model.status() == LP::Status::ABORTED_NO_FEASIBLE)
	    fvs = heurFvs;
	//std::cerr << "bound = " << model.bestBound()
	//	  << " lb = " << lb << " fvs.size() = " << fvs.size() << std::endl;
	assert(lb < fvs.size());
	break;
        }
    default:
	abort();
    }

    delete &model;
    return {lb, fvs};
}

int main(int argc, char* argv[]) {
    double tStart = now();
    bool statsOnly = false;
    int opt;
    while ((opt = getopt(argc, argv, "st:v")) != -1) {
        switch (opt) {
        case 's':
            statsOnly = true;
            break;
        case 't':
            timeLimit = std::stod(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        default:
            usage(argv[0]);
            exit(1);
            break;
        }
    }

    Graph g0;
    std::map<Vertex, std::string> vertexNames;
    std::tie(g0, vertexNames) = Graph::read(std::cin);

    size_t best_lb = 0;
    size_t best_ub = g0.numVertices();

    Graph g = g0;
    std::vector<Vertex> fvs;

    //for (Edge e : g.edges()) std::cout << vertexNames[e.u] << ' ' << vertexNames[e.v] << std::endl;
    //for (Edge e : g.edges()) std::cout << e.u << ' ' << e.v << std::endl;
    std::vector<Graph> reducedComponents;
    for (Graph gc : g.connectedComponents()) {
	Graph rgc;
	std::vector<Vertex> fvsc;
	std::tie(rgc, fvsc) = reduce(gc, vertexNames);
	for (Vertex u : fvsc)
	    fvs.push_back(u);
	for (Graph gcc : rgc.connectedComponents())
	    reducedComponents.push_back(gcc);
    }

    size_t selfLoops = fvs.size();

    std::vector<size_t> lbs;
    std::vector<std::vector<Vertex>> fvss;
    if (reducedComponents.empty()) {
	std::cerr << "reduced-n: " << 0 << std::endl
		  << "reduced-m: " << 0 << std::endl
		  << "reduced-k-lb: " << 0 << std::endl
		  << "reduced-k-ub: " << 0 << std::endl;
    } else {
	std::sort(reducedComponents.begin(), reducedComponents.end(),
		  [](const Graph& g1, const Graph& g2) {
		      return g1.numVertices() < g2.numVertices();
		  });
	std::cerr << "reduced-n: " << reducedComponents.back().numVertices() << std::endl;
	std::cerr << "reduced-m: " << reducedComponents.back().numEdges() << std::endl;

#if 0
	for (const auto& g : reducedComponents) {
	    std::cerr << "component:\n";
	    for (Edge e : g.edges()) std::cerr << e.u << ' ' << e.v << std::endl;
	    std::cerr << "cycles:\n";
	    for (auto cycle : chordlessCycles(g)) {
		bool first = true;
		for (auto u : cycle) {
		    if (!first)
			std::cerr << ' ';
		    std::cerr << u;
		    first = false;
		}
		std::cout << std::endl;
	    }
	}
	//exit(0);
#endif

	for (const auto& gc : reducedComponents)
	    lbs.push_back(fvsLB(gc));
	std::cerr << "reduced-k-lb: " << lbs.back() << std::endl;

	for (const auto& gc : reducedComponents)
	    fvss.push_back(fvsHeuristic(gc));
	std::cerr << "reduced-k-ub: " << fvss.back().size() << std::endl;
    }

    for (size_t i = 0; i < reducedComponents.size(); ++i) {
	size_t lb = selfLoops;
	for (auto b : lbs)
	    lb += b;
	size_t ub = selfLoops;
	for (const auto& f : fvss)
	    ub += f.size();
	if (lb > best_lb) {
	    std::cerr << "lb: " << lb << std::endl;
	    best_lb = lb;
	}
	if (ub < best_ub) {
	    std::cerr << "ub: " << ub << std::endl;
	    best_ub = ub;
	}

	size_t other_lb = lb - lbs[i];
	size_t other_ub = ub - fvss[i].size();
	size_t best_comp_lb = i == reducedComponents.size() - 1 ? lbs.back() : 1e10;
	size_t best_comp_ub = i == reducedComponents.size() - 1 ? fvss.back().size() : 0;

	const auto& gc = reducedComponents[i];
	const auto& fvscHeuristic = fvss[i];
	size_t lbc;
	std::vector<Vertex> fvscILP;
#if 1
	std::tie(lbc, fvscILP) = fvsHS(gc, fvscHeuristic,
				       other_lb, other_ub, best_lb, best_ub,
				       best_comp_lb, best_comp_ub);
#else
	std::tie(lbc, fvscILP) = fvsILP(gc, fvscHeuristic,
					other_lb, other_ub, best_lb, best_ub,
					best_comp_lb, best_comp_ub);
#endif
	if (verbose)
	    fprintf(stderr, "%6zd %7zd %5zd %5zd %5zd\n",
		    gc.numVertices(), gc.numEdges(), lbc, fvscILP.size(), fvscHeuristic.size());
	auto fvsc = fvscILP.size() < fvscHeuristic.size() ? fvscILP : fvscHeuristic;
	if (i == reducedComponents.size() - 1)
	    std::cerr << "reduced-k: " << fvsc.size() << std::endl;
	for (Vertex u : fvsc)
	    fvs.push_back(u);
	lbs[i] = lbc;
	fvss[i] = fvsc;
    }
    double tStop = now();

    Graph g2 = g0;
    for (Vertex u : fvs)
	g2.deleteVertex(u);
    auto cycles = findCycles(g2);
    if (!cycles.empty()) {
	std::set<std::set<Vertex>> done;
	std::cerr << "internal error: cycles remain\n";
	std::cerr << "fvs:";
	for (Vertex u : fvs)
	    std::cerr << ' ' << vertexNames[u];
	std::cerr << std::endl;
	for (auto cycle : cycles) {
	    std::set<Vertex> cs(cycle.begin(), cycle.end());
	    if (done.find(cs) != done.end())
		continue;
	    std::cerr << "cycle:";
	    for (Vertex u : cs)
		std::cerr << ' ' << vertexNames[u];
	    std::cerr << std::endl;
	    done.insert(cs);
	}
	exit(1);
    }

    size_t lb = selfLoops;
    for (auto b : lbs)
	lb += b;

    if (statsOnly) {
	printf("%6zd %7zd "
	       "%5zd %5zd "
	       "%6zd %7zd "
	       "%5zd %5zd "
	       "%9.2f\n",
	       g0.numVertices(), g0.numEdges(),
	       lb, fvs.size(),
	       reducedComponents.back().numVertices(), reducedComponents.back().numEdges(),
	       lbs.back(), fvss.back().size(),
	       tStop - tStart);
    } else {
	for (Vertex u : fvs)
	    std::cout << vertexNames[u] << std::endl;
    }

    return 0;
}
