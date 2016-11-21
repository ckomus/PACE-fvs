/* commode -- find significant protein complex graph models
   Copyright (C) 2014  Falk Hüffner

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

#include <stdexcept>
#include <vector>

#include <cmath>
#include <assert.h>

#include "LP.hh"
#ifdef HAVE_CPLEX
# include "LPCplex.hh"
#endif
#ifdef HAVE_GUROBI
# include "LPGurobi.hh"
#endif

namespace LP {

Model* newModel() {
    return newModel({Solver::GUROBI, Solver::CPLEX});
}

Model* newModel(Solver solver) {
    (void) solver;
#ifdef HAVE_CPLEX
    if (solver == Solver::CPLEX)
	return new CplexModel();
#endif
#ifdef HAVE_GUROBI
    if (solver == Solver::GUROBI)
	return new GurobiModel();
#endif
    return nullptr;
}

Model* newModel(const std::vector<Solver>& solvers) {
    for (Solver solver : solvers) {
	auto m = newModel(solver);
	if (m)
	    return m;
    }
    return nullptr;
}

double Expr::value() const {
    double x = constant();
    for (auto& p : coeff_)
	x += Var(*m_, p.first).value() * p.second;
    return x;
}

Expr Expr::operator-() const {
    Expr e = *this;
    for (auto& p : e.coeff_)
	p.second = -p.second;
    e.constant_ = -e.constant_;
    return e;
}

Expr& Expr::operator+=(const Expr& e) {
    if (isNull())
	return *this = e;
    if (e.isNull())
	return *this;
    if (&e.model() != &model())
	throw std::logic_error("cannot add expressions from different models");
    for (auto& p : e.coeff_) {
	auto q = coeff_.find(p.first);
	if (q == coeff_.end()) {
	    coeff_[p.first] = p.second;
	} else {
	    double c = q->second + p.second;
	    if (c == 0.0)
		coeff_.erase(q);
	    else
		q->second = c;
	}
    }
    constant_ += e.constant_;
    return *this;
}

Expr& Expr::operator-=(const Expr& e) {
    if (isNull())
	return *this = -e;
    if (e.isNull())
	return *this;
    if (&e.model() != &model())
	throw std::logic_error("cannot add expressions from different models");
    for (auto& p : e.coeff_) {
	auto q = coeff_.find(p.first);
	if (q == coeff_.end()) {
	    coeff_[p.first] = -p.second;
	} else {
	    double c = q->second - p.second;
	    if (c == 0.0)
		coeff_.erase(q);
	    else
		q->second = c;
	}
    }
    constant_ -= e.constant_;
    return *this;
}

Expr& Expr::operator*=(double c) {
    if (c == 0) {
	coeff_ = {};
    } else {
	for (auto& p : coeff_)
	    p.second *= c;
    }
    constant_ *= c;
    return *this;
}

std::ostream& operator<<(std::ostream& out, const LP::Expr& expr) {
    bool start = true;
    for (auto term : expr) {
	double c = term.coefficient();
	assert(c != 0);
	if (start) {
	    if (c < 0)
		out << "- ";
	} else {
	    out << (c < 0 ? " - " : " + ");
	}
	if (c != 1 && c != -1)
	    out << std::abs(c) << ' ';
	out << term.var().name();
	start = false;
    }
    if (expr.constant() != 0) {
	if (start) {
	    if (expr.constant() < 0)
		out << "- ";
	} else {
	    out << (expr.constant() < 0 ? " - " : " + ");
	}
	out << std::abs(expr.constant());
    } else if (expr.numCoefficients() == 0) {
	out << '0';
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const LP::ConstrExpr& constr) {
    Expr e = constr.lhs();
    double c = e.constant();
    e -= c;
    out << e;
    switch (constr.sense()) {
    case ConstrExpr::Sense::LE: out << " <= "; break;
    case ConstrExpr::Sense::EQ: out << " = ";  break;
    case ConstrExpr::Sense::GE: out << " >= "; break;
    }
    return out << -c;
}

std::ostream& operator<<(std::ostream& out, const LP::Model& model) {
    out << (model.objectiveSense() == Model::Sense::MINIMIZE ? "minimize\n" : "maximize\n")
	<< model.objective()
	<< "\nsubject to\n";
    for (auto c = model.constraintsBegin(); c != model.constraintsEnd(); ++c)
	out << (*c).constrExpr() << '\n';
    out << "bounds\n";
    for (auto i = model.varsBegin(); i != model.varsEnd(); ++i) {
	Var v = *i;
	out << v.lb() << " <= " << v.name() << " <= " << v.ub() << '\n';
    }
    out << "binaries\n";
    for (auto i = model.varsBegin(); i != model.varsEnd(); ++i) {
	Var v = *i;
	if (v.isBool())
	    out << v.name() << '\n';
    }
    out << "generals\n";
    for (auto i = model.varsBegin(); i != model.varsEnd(); ++i) {
	Var v = *i;
	if (v.isIntegral() && !v.isBool())
	    out << v.name() << '\n';
    }
    return out << "end" << std::endl;
}

}
