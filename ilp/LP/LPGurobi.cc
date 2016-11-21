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

#include <cassert>
#include <cmath>

#include <functional>
#include <string>
#include <vector>

#include "LPGurobi.hh"

namespace LP {

#define CHECKSTATUS(fun, args) do {					\
	if (fun args)							\
	    error(std::string(__func__) + ": " #fun);			\
    } while (0)

GurobiModel::GurobiModel(Sense sense) : addedVars_(0), addedConstrs_(0) {
    if (!env_)
	if (GRBloadenv(&env_, nullptr))
	    error("GurobiModel::GurobiModel: GRBloadenv");
    if (GRBnewmodel(env_, &lp_, nullptr, 0, nullptr, nullptr, nullptr, nullptr, nullptr))
	error("GurobiModel::GurobiModel: GRBnewmodel");
    if (sense == Sense::MAXIMIZE)
	setObjectiveSense(Sense::MAXIMIZE);
    setVerbose(false);
#ifdef LP_SINGLE_THREAD
    GRBsetintparam(env_, "Threads", 1);
#endif
}

GurobiModel::~GurobiModel() {
    if (GRBfreemodel(lp_))
	error("GurobiModel::~GurobiModel: GRBfreemodel");
}

void GurobiModel::setVerbose(bool verbose) {
    if (GRBsetintparam(env(), "OutputFlag", verbose))
	error("GurobiModel::setVerbose: GRBsetintparam(OutputFlag)");
}

// querying

size_t GurobiModel::numVariables() const {
    int numVars;
    if (GRBgetintattr(lp_, "NumVars", &numVars))
	error("GurobiModel::numVariables: GRBgetintattr(NumVars)");
    return numVars;
}

size_t GurobiModel::numConstraints() const {
    int numConstrs;
    if (GRBgetintattr(lp_, "NumConstrs", &numConstrs))
	error("GurobiModel::numConstraints: GRBgetintattr(NumConstrs)");
    return numConstrs;
}

Expr GurobiModel::objective() const {
    update();
    size_t n = numVariables();
    double* obj = new double[n];
    if (GRBgetdblattrarray(lp_, "Obj", 0, n, obj))
	error("GurobiModel::objective: GRBgetdblattrarray(Obj)");

    Expr e;
    for (size_t i = 0; i < n; i++)
	if (obj[i])
	    e += obj[i] * makeVar(i);
    delete[] obj;

    double objoffset;
    if (GRBgetdblattr(lp_, "ObjCon", &objoffset))
	error("GurobiModel::objective: GRBgetdblattr(ObjCon)");
    e += objoffset;

    return e;
}

Model::Sense GurobiModel::objectiveSense() const {
    int sense;
    if (GRBgetintattr(lp_, "ModelSense", &sense))
	error("GurobiModel::objectiveSense: GRBgetintattr(ModelSense)");

    switch (sense) {
    case +1: return Sense::MINIMIZE; break;
    case -1: return Sense::MAXIMIZE; break;
    }
    error("GurobiModel::objectiveSense: unknown objective sense");
    return Sense::MINIMIZE; // suppress warning
}

double GurobiModel::objectiveValue() const {
    double objval;
    if (GRBgetdblattr(lp_, "ObjVal", &objval))
	error("GurobiModel::varName: GRBgetdblattr(ObjVal)");
    return objval;
}

double GurobiModel::bestBound() const {
    double bestobjval;
    CHECKSTATUS(GRBgetdblattr, (lp_, "ObjBound", &bestobjval));
    return bestobjval;
}

std::string GurobiModel::varName(int v) const {
    char *name;
    if (GRBgetstrattrelement(lp_, "VarName", v, &name))
	error("GurobiModel::varName: GRBgetstrattrelement(VarName)");
    return name;
}

ConstrExpr GurobiModel::constrExpr(const Constraint& c) const {
    int row = c.id();
    int nzcnt;
    int rmatbeg[1];
    if (GRBgetconstrs(lp_,  &nzcnt, nullptr, nullptr, nullptr, row, 1))
	error("GurobiModel::constrExpr: GRBgetconstrs");
    int* rmatind = new int[nzcnt];
    double* rmatval = new double[nzcnt];
    if (GRBgetconstrs(lp_, &nzcnt, rmatbeg, rmatind, rmatval, row, 1))
	error("GurobiModel::constrExpr: GRBgetconstrs");
    assert(rmatbeg[0] == 0);

    Expr e;
    for (size_t i = 0; i < size_t(nzcnt); i++)
	e += rmatval[i] * makeVar(rmatind[i]);
    delete[] rmatind;
    delete[] rmatval;

    double rhs;
    if (GRBgetdblattrelement(lp_, "RHS", row, &rhs))
	error("GurobiModel::constrExpr: GRBgetdblattrelement(RHS)");
    e -= rhs;

    char gurobiSense;
    if (GRBgetcharattrelement(lp_, "Sense", row, &gurobiSense))
	error("GurobiModel::constrExpr: GRBgetcharattrelement(Sense)");
    ConstrExpr::Sense sense;
    switch (gurobiSense) {
    case '<': sense = ConstrExpr::Sense::LE; break;
    case '=': sense = ConstrExpr::Sense::EQ; break;
    case '>': sense = ConstrExpr::Sense::GE; break;
    default: error("GurobiModel::constrExpr: unknown column type"); break;
    }
    return ConstrExpr(e, sense);
}

// modifying

Var GurobiModel::doAddVar(double lb, double ub, double obj, const std::string& name, char type) {
    int v = numVariables() + addedVars_++;
    if (GRBaddvar(lp_, 0, nullptr, nullptr, obj, lb, ub, type, const_cast<char*>(name.c_str())))
	error("GurobiModel::doAddVar: GRBaddvar");
    return makeVar(v);
}

Var GurobiModel::addVar(double lb, double ub, double obj, const std::string& name) {
    return doAddVar(lb, ub, obj, name,  GRB_CONTINUOUS);
}

Var GurobiModel::addIntVar(double lb, double ub, double obj, const std::string& name) {
    return doAddVar(lb, ub, obj, name, lb == 0.0 && ub == 1.0 ? GRB_BINARY : GRB_INTEGER);
}

static char gurobiSense(ConstrExpr::Sense sense) {
    switch (sense) {
    case ConstrExpr::Sense::LE: return GRB_LESS_EQUAL;
    case ConstrExpr::Sense::EQ: return GRB_EQUAL;
    case ConstrExpr::Sense::GE: return GRB_GREATER_EQUAL;
    }
    abort();
}

Constraint GurobiModel::addConstraint(const ConstrExpr& c) {
    if (addedVars_)
	update();
    int row = numConstraints() + addedConstrs_++;
    std::vector<int> rmatind;
    std::vector<double> rmatval;
    for (auto term : c.lhs()) {
	rmatind.push_back(term.var().id());
	rmatval.push_back(term.coefficient());
    }
    int nzcnt = rmatind.size();
    double rhs = -c.lhs().constant();
    if (GRBaddconstr(lp_, nzcnt, &rmatind[0], &rmatval[0], gurobiSense(c.sense()), rhs, nullptr))
	error("GurobiModel::addConstraint: GRBaddconstr");
    return makeConstraint(row);
}

Constraint GurobiModel::addConditionalConstraint(const Expr& condition, const ConstrExpr& constraint) {
    condition.extractBoolVar();
    double lhsLB = constraint.lhs().constant();
    double lhsUB = constraint.lhs().constant();
    for (auto term : constraint.lhs()) {
	if (term.coefficient() >= 0) {
	    lhsLB += term.coefficient() * term.var().lb();
	    lhsUB += term.coefficient() * term.var().ub();
	} else {
	    lhsLB += term.coefficient() * term.var().ub();
	    lhsUB += term.coefficient() * term.var().lb();
	}
    }

    // FIXME: for equality, we only return one of the resulting constraints.
    if (constraint.sense() == ConstrExpr::Sense::LE || constraint.sense() == ConstrExpr::Sense::EQ) {
	Constraint result = addConstraint(constraint.lhs() <= 0 + lhsUB * !condition);
	if (constraint.sense() == ConstrExpr::Sense::LE)
	    return result;
    }
    // GE or EQ
    return addConstraint(constraint.lhs() >= 0 + lhsLB * !condition);
}

void GurobiModel::solve() {
    if (GRBoptimize(lp_))
	error("GurobiModel::solve: GRBoptimize");
}

Status GurobiModel::status() const {
    int state;
    if (GRBgetintattr(lp_, "Status", &state))
	error("GurobiModel::status: GRBgetintattr(Status)");
    switch (state) {
    case GRB_OPTIMAL:     return Status::OPTIMAL;                 break;
    case GRB_INFEASIBLE:  return Status::INFEASIBLE;		  break;
    case GRB_INF_OR_UNBD: return Status::INFEASIBLE_OR_UNBOUNDED; break;
    case GRB_INTERRUPTED: case GRB_TIME_LIMIT: {
	int solCount;
	if (GRBgetintattr(lp_, "SolCount", &solCount))
	    error("GurobiModel::status: GRBgetintattr(SolCount)");
	if (solCount)
	    return Status::ABORTED_FEASIBLE;
	else
	    return Status::ABORTED_NO_FEASIBLE;
	break;
    }
    default:
	throw std::runtime_error("GurobiModel::status: unknown Gurobi status " + std::to_string(state));
	break;
    }
}

bool GurobiModel::varIsIntegral(const Var& v) const {
    update();
    char type = varCharAttr(v, "VType");
    switch (type) {
    case GRB_INTEGER: case GRB_BINARY: return true;  break;
    case GRB_CONTINUOUS:               return false; break;
    default:
	throw std::runtime_error("GurobiModel::status: unknown Gurobi variable type " + type);
	break;
    }
}

double GurobiModel::varlb(const Var& v) const {
    return varDoubleAttr(v, "LB");
}

double GurobiModel::varub(const Var& v) const {
    return varDoubleAttr(v, "UB");
}

double GurobiModel::varValue(const Var& v) const {
    double x = varDoubleAttr(v, "X");
    // Gurobi sometimes returns slightly nonintegral values even for
    // integer variables. To avoid messy logic, we properly round them
    // here.
    char type = varCharAttr(v, "VType");
    switch (type) {
    case GRB_INTEGER: case GRB_BINARY:
	x = std::round(x);
	break;
    case GRB_CONTINUOUS:
	break;
    default:
	throw std::runtime_error("GurobiModel::varValue: unknown Gurobi variable type " + type);
	break;
    }
    return x;
}

double GurobiModel::dualValue(const Constraint& c) const {
    double x;
    if (GRBgetdblattrelement(lp_, "Pi", c.id(), &x))
	error("GurobiModel::dualValue: GRBgetdblattrelement(RC)");
    return x;
}

void GurobiModel::varChangeBounds(const Var& v, double lb, double ub) const {
    if (GRBsetdblattrelement(lp_, "LB", v.id(), lb))
	error("GurobiModel::varChangeBounds: GRBsetdblattrelement(LB)");
    if (GRBsetdblattrelement(lp_, "UB", v.id(), ub))
	error("GurobiModel::varChangeBounds: GRBsetdblattrelement(UB)");
}

void GurobiModel::varSetIntegral(const Var& v, bool integral) const {
    if (GRBsetcharattrelement(lp_, "VType", v.id(), integral ? 'I' : 'B'))
	error("GurobiModel::varSetIntegral: GRBsetcharattrelement(VType)");
}

void GurobiModel::setCoefficient(const Constraint& c, const Var& v, double x) const {
    if (addedVars_ || addedConstrs_)
	update();
    int cind[] = { c.id() };
    int vind[] = { v.id() };
    double val[] = { x };
    if (GRBchgcoeffs(lp_, 1, cind, vind, val))
	error("GurobiModel::setCoefficient: GRBchgcoeffs");
}

void GurobiModel::setObjective(const Expr& obj) {
    update();
    size_t n = numVariables();
    std::vector<double> coeff(n);
    for (auto term : obj)
	coeff[term.var().id()] = term.coefficient();
    if (GRBsetdblattrarray(lp_, "Obj", 0, n, &coeff[0]))
	error("GurobiModel::setObjective::GRBsetdblattrarray");
    setObjectiveConstant(obj.constant());
}

void GurobiModel::setObjectiveCoeff(const Var& v, double c) {
    if (GRBsetdblattrelement(lp_, "Obj", v.id(), c))
	error("GurobiModel::setObjectiveCoeff: GRBsetdblattrelement(Obj)");
}

void GurobiModel::setObjectiveSense(Sense s) {
    if (GRBsetintattr(lp_, "ModelSense", s == Sense::MINIMIZE ? +1 : -1))
	error("GurobiModel::setObjectiveSense: GRBsetintattrelement(ModelSense)");
}

void GurobiModel::setObjectiveConstant(double x) {
    if (GRBsetdblattr(lp_, "ObjCon", x))
	error("GurobiModel::setObjectiveConstant::GRBsetdblattr");
}

void GurobiModel::MIPstart(const std::map<Var, double>& vars) {
    std::vector<int> indices;
    std::vector<double> values;
    for (auto p : vars) {
        indices.push_back(p.first.id());
        values.push_back(p.second);
    }
    CHECKSTATUS(GRBsetdblattrlist, (lp_, "Start", indices.size(), indices.data(), values.data()));
}

void GurobiModel::setEmphasis(Emphasis e) {
    int emph = GRB_MIPFOCUS_BALANCED;
    switch (e) {
    case Emphasis::BALANCED:         emph = GRB_MIPFOCUS_BALANCED;    break;
    case Emphasis::FIND_FEASIBLE:    emph = GRB_MIPFOCUS_FEASIBILITY; break;
    case Emphasis::PROVE_OPTIMALITY: emph = GRB_MIPFOCUS_OPTIMALITY;  break;
    case Emphasis::MOVE_BOUND:       emph = GRB_MIPFOCUS_BESTBOUND;   break;
    }
    if (GRBsetintparam(env_, GRB_INT_PAR_MIPFOCUS, emph))
	error("GurobiModel::setEmphasis");
}

class GurobiCallbackSolver : public CallbackSolver {
public:
    GurobiCallbackSolver(GRBmodel* lp, void* cbdata, int where, bool isIntegral) :
	lp_(lp), cbdata_(cbdata), where_(where), isIntegral_(isIntegral) {
	int numColumns;
	CHECKSTATUS(GRBgetintattr, (lp_, "NumVars", &numColumns));
	x_.resize(numColumns);
	CHECKSTATUS(GRBcbget, (cbdata, where, isIntegral ? GRB_CB_MIPSOL_SOL : GRB_CB_MIPNODE_REL, &x_[0]));
	if (isIntegral) {
	    std::vector<char> colTypes(x_.size());
	    CHECKSTATUS(GRBgetcharattrarray, (lp_, "VType", 0, colTypes.size(), &colTypes[0]));
	    for (std::size_t i = 0; i < x_.size(); i++)
		if (colTypes[i] != GRB_CONTINUOUS)
		    x_[i] = intRound(x_[i]);
	}
    }
    double value(const Var& v) const { return x_[v.id()]; }
    double bestBound() const {
	int what;
	switch(where_) {
	case GRB_CB_MIP: what = GRB_CB_MIP_OBJBND; break;
	case GRB_CB_MIPSOL: what = GRB_CB_MIPSOL_OBJBND; break;
	case GRB_CB_MIPNODE: what = GRB_CB_MIPNODE_OBJBND; break;
	default: {
	    int sense;
	    CHECKSTATUS(GRBgetintattr, (lp_, "ModelSense", &sense));
	    if (sense > 0)
		return -LP::inf;
	    else
		return +LP::inf;
	}
	}
	double objbnd;
	CHECKSTATUS(GRBcbget, (cbdata_, where_, what, &objbnd));
	//fprintf(stderr, "best bound: %f\n", objbnd);
	if (objbnd >= GRB_INFINITY)
	    objbnd = +LP::inf;
	if (objbnd <= -GRB_INFINITY)
	    objbnd = -LP::inf;
	return objbnd;
    }
    double bestObjectiveValue() const {
	int what;
	switch(where_) {
	case GRB_CB_MIP: what = GRB_CB_MIP_OBJBST; break;
	case GRB_CB_MIPSOL: what = GRB_CB_MIPSOL_OBJBST; break;
	case GRB_CB_MIPNODE: what = GRB_CB_MIPNODE_OBJBST; break;
	default: {
	    int sense;
	    CHECKSTATUS(GRBgetintattr, (lp_, "ModelSense", &sense));
	    if (sense > 0)
		return +LP::inf;
	    else
		return -LP::inf;
	}
	}
	double objbst;
	CHECKSTATUS(GRBcbget, (cbdata_, where_, what, &objbst));
	return objbst;
    }
    void addLazyConstraint(const ConstrExpr& constr) { addConstraint(constr, GRBcblazy); }
    void addUserCut(const ConstrExpr& constr) { addConstraint(constr, GRBcbcut); }

private:
    void addConstraint(const ConstrExpr& constr,
		       std::function<int(void*, int, const int*, const double*, char, double)> f) {
	std::vector<int> indices;
	std::vector<double> values;
	for (auto term : constr.lhs()) {
	    indices.push_back(term.var().id());
	    values.push_back(term.coefficient());
	}
	double rhs = -constr.lhs().constant();
	CHECKSTATUS(f, (cbdata_, indices.size(), &indices[0], &values[0], gurobiSense(constr.sense()), rhs));
    }

    void error(const std::string& msg) const {
	fprintf(stderr, "|||%s|||\n", GRBgeterrormsg(GRBgetenv(lp_)));
	throw std::runtime_error(msg + ": " + GRBgeterrormsg(GRBgetenv(lp_)) + "!!!");
    }
    GRBmodel* lp_;
    void* cbdata_;
    int where_;
    bool isIntegral_;
    std::vector<double> x_;
};


static int gurobiCallback(GRBmodel* lp, void* cbdata, int where, void* usrdata) {
    GurobiModel* model = reinterpret_cast<GurobiModel*>(usrdata);
    if (where == GRB_CB_MIPSOL && (model->callbackType() & LAZY_CONSTRAINTS)) {
	GurobiCallbackSolver callbackSolver(lp, cbdata, where, true);
	CallbackModel callbackModel(callbackSolver);
	model->callback(callbackModel);
    } else if (where == GRB_CB_MIPNODE && (model->callbackType() & USER_CUTS)) {
	int status;
	if (GRBcbget(cbdata, where, GRB_CB_MIPNODE_STATUS, &status))
            throw std::runtime_error(std::string(__func__) + ": GRBcbget(GRB_CB_MIPNODE_STATUS)");
        if (status == GRB_LOADED || status == GRB_INFEASIBLE || status == GRB_CUTOFF)
            return 0;
        if (status != GRB_OPTIMAL)
	    throw std::runtime_error(std::string(__func__) + ": unknown status");
	GurobiCallbackSolver callbackSolver(lp, cbdata, where, false);
	CallbackModel callbackModel(callbackSolver);
	model->callback(callbackModel);
    }
    return 0;
}

void GurobiModel::setCallback(Callback& callback, CallbackType type) {
    callback_ = &callback;
    callbackType_ = type;
    if (type & LAZY_CONSTRAINTS)
	CHECKSTATUS(GRBsetintparam, (env(), "LazyConstraints", 1));
    if (type & USER_CUTS)
	CHECKSTATUS(GRBsetintparam, (env(), "PreCrush", 1));
    CHECKSTATUS(GRBsetcallbackfunc, (lp_, gurobiCallback, this));
}

void GurobiModel::setTimeLimit(double seconds) {
    CHECKSTATUS(GRBsetdblparam, (env(), "TimeLimit", seconds));
}

int GurobiModel::varIntAttr(const Var& v, const char* attr) const {
    int x;
    if (GRBgetintattrelement(lp_, attr, v.id(), &x))
	error("GRBgetintattr(" + std::string(attr) + ")");
    return x;
}

char GurobiModel::varCharAttr(const Var& v, const char* attr) const {
    char x;
    if (GRBgetcharattrelement(lp_, attr, v.id(), &x))
	error("GRBgetcharattr(" + std::string(attr) + ")");
    return x;
}

double GurobiModel::varDoubleAttr(const Var& v, const char* attr) const {
    double x;
    if (GRBgetdblattrelement(lp_, attr, v.id(), &x))
	error("GRBgetdblattr(" + std::string(attr) + ")");
    return x;
}

const char* GurobiModel::varStringAttr(const Var& v, const char* attr) const {
    char *x;
    if (GRBgetstrattrelement(lp_, attr, v.id(), &x))
	error("GRBgetdblattr(" + std::string(attr) + ")");
    return x;
}

void GurobiModel::error(const std::string& msg) const {
    throw std::runtime_error(msg + ": " + GRBgeterrormsg(env()) + "!!!");
}

GRBenv* GurobiModel::env_;

}
