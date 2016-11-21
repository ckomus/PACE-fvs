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
#include <stdexcept>
#include <vector>
#include <cmath>

#include <fcntl.h>
#include <unistd.h>

#include "LPCplex.hh"

namespace LP {

#define CHECKSTATUS(fun, args) do {					\
	int CHECKSTATUS_status = fun args;				\
	if (CHECKSTATUS_status)						\
	    error(CHECKSTATUS_status, std::string(__func__) + ": " #fun); \
    } while (0)


CplexModel::CplexModel(Sense sense) : objConstant_(0) {
    // suppress annoying "license manager" output
    int backup, redirection;
    fflush(stderr);
    backup = dup(2);
    redirection = open("/dev/null", O_WRONLY);
    dup2(redirection, 2);
    close(redirection);
    int status;
    env_ = CPXopenCPLEX(&status);
    fflush(stderr);
    dup2(backup, 2);
    close(backup);
    if (!env_)
	error(status, "CplexModel::CplexModel");

    lp_ = CPXcreateprob(env_, &status, "model");
    if (!lp_)
	error(status, "CplexModel::CplexModel");

    if (sense == Sense::MAXIMIZE)
	CPXchgobjsen(env_, lp_, CPX_MAX);
#ifdef LP_SINGLE_THREAD
    CPXsetintparam(env_, CPXPARAM_Threads, 1);
#endif
}

CplexModel::~CplexModel() {
    int status = CPXcloseCPLEX(&env_);
    if (status)
	error(status, "CplexModel::~CplexModel");
}

void CplexModel::setVerbose(bool verbose) {
    int status;
    if ((status = CPXsetintparam(env_, CPXPARAM_ScreenOutput, verbose ? CPX_ON : CPX_OFF)))
	error(status, "CplexModel::setVerbose: CPXsetintparam(CPXPARAM_ScreenOutput)");
}

// querying

size_t CplexModel::numVariables() const {
    return CPXgetnumcols(env_, lp_);
}

size_t CplexModel::numConstraints() const {
    return CPXgetnumrows(env_, lp_);
}

Expr CplexModel::objective() const {
    size_t n = numVariables();
    Expr e;
    int status;
    if (n) {
	double* obj = new double[n];
	status = CPXgetobj(env_, lp_, obj, 0, n - 1);
	if (status)
	    error(status, "CplexModel::objective: CPXgetobj");
	for (size_t i = 0; i < n; i++)
	    if (obj[i])
		e += obj[i] * makeVar(i);
	delete[] obj;
    }
    double objoffset;
    if ((status = CPXgetobjoffset(env_, lp_, &objoffset)))
	error(status, "CplexModel::objective: CPXgetobjoffset");
    e += objoffset;

    return e;
}

Model::Sense CplexModel::objectiveSense() const {
    switch (CPXgetobjsen(env_, lp_)) {
    case CPX_MIN: return Sense::MINIMIZE; break;
    case CPX_MAX: return Sense::MAXIMIZE; break;
    }
    throw std::runtime_error("CplexModel::objectiveSense: unknown objective sense");
}

double CplexModel::objectiveValue() const {
    double objval;
    int status = CPXgetobjval(env_, lp_, &objval);
    if (status)
	error(status, "CplexModel::objectiveValue");
    return objval + objConstant_;
}

double CplexModel::bestBound() const {
    // FIXME: round up if integral objective
    double bestobjval;
    CHECKSTATUS(CPXgetbestobjval, (env_, lp_, &bestobjval));
    return bestobjval + objConstant_;
}

std::string CplexModel::varName(int v) const {
    char* name[1];
    char namestore[4096];
    int surplus;
    int status = CPXgetcolname(env_, lp_, name, namestore, sizeof namestore, &surplus, v, v);
    if (status == CPXERR_NEGATIVE_SURPLUS) {
	size_t sizeof_namestore2 = sizeof namestore + (-surplus);
	char* namestore2 = new char[sizeof_namestore2];
	status = CPXgetcolname(env_, lp_, name, namestore2, sizeof_namestore2, &surplus, v, v);
	delete[] namestore2;
    }
    if (status)
	error(status, "CplexModel::varName");
    return name[0];
}

ConstrExpr CplexModel::constrExpr(const Constraint& c) const {
    int row = c.id();
    int nzcnt;
    int rmatbeg[1];
    int surplus;
    int status = CPXgetrows(env_, lp_, &nzcnt, rmatbeg, nullptr, nullptr, 0, &surplus, row, row);
    //if (status != CPXERR_NEGATIVE_SURPLUS)
    //    error(status, "CplexModel::constrExpr");
    size_t n = -surplus;
    int* rmatind = new int[n];
    double* rmatval = new double[n];
    if ((status = CPXgetrows(env_, lp_, &nzcnt, rmatbeg, rmatind, rmatval, n, &surplus, row, row)))
	error(status, "CplexModel::constrExpr: CPXgetrows");
    assert(surplus == 0);
    assert(rmatbeg[0] == 0);

    Expr e;
    for (size_t i = 0; i < n; i++)
	e += rmatval[i] * makeVar(rmatind[i]);
    delete[] rmatind;
    delete[] rmatval;

    double rhs;
    if ((status = CPXgetrhs(env_, lp_, &rhs, row, row)))
	error(status, "CplexModel::constrExpr: CPXgetsense");
    e -= rhs;

    char cplex_sense;
    if ((status = CPXgetsense(env_, lp_, &cplex_sense, row, row)))
	error(status, "CplexModel::constrExpr: CPXgetsense");
    ConstrExpr::Sense sense;
    switch (cplex_sense) {
    case 'L': sense = ConstrExpr::Sense::LE; break;
    case 'E': sense = ConstrExpr::Sense::EQ; break;
    case 'G': sense = ConstrExpr::Sense::GE; break;
    default: error(0, "CplexModel::constrExpr: unknown column type"); break;
    }
    return ConstrExpr(e, sense);
}

// modifying

Var CplexModel::doAddVar(double lb, double ub, double obj, const std::string& name, const char* type) {
    int v =  CPXgetnumcols(env_, lp_);
    char* names[1]{const_cast<char*>(name.c_str())};
    int status = CPXnewcols(env_, lp_, 1, &obj, &lb, &ub, type, names);
    if (status)
	error(status, "CplexModel::addVar");
    return makeVar(v);
}

Var CplexModel::addVar(double lb, double ub, double obj, const std::string& name) {
    return doAddVar(lb, ub, obj, name, nullptr);
}

Var CplexModel::addIntVar(double lb, double ub, double obj, const std::string& name) {
    return doAddVar(lb, ub, obj, name, lb == 0.0 && ub == 1.0 ? "B" : "I");
}

static char cplexConstrSense(ConstrExpr::Sense sense) {
    switch (sense) {
    case ConstrExpr::Sense::LE: return 'L';
    case ConstrExpr::Sense::EQ: return 'E';
    case ConstrExpr::Sense::GE: return 'G';
    }
    abort();
}

Constraint CplexModel::addConstraint(const ConstrExpr& c) {
    int row = CPXgetnumrows(env_, lp_);
    std::vector<int> rmatind;
    std::vector<double> rmatval;
    for (auto term : c.lhs()) {
	rmatind.push_back(term.var().id());
	rmatval.push_back(term.coefficient());
    }
    int nzcnt = rmatind.size();
    double rhs = -c.lhs().constant();
    char sense = cplexConstrSense(c.sense());
    int rmatbeg[1] = {0};
    int status = CPXaddrows(env_, lp_, 0, 1, nzcnt, &rhs, &sense, rmatbeg,
			    &rmatind[0], &rmatval[0], nullptr, nullptr);
    if (status)
	error(status, "CplexModel::addConstraint");
    return makeConstraint(row);
}

Constraint CplexModel::addConditionalConstraint(const Expr& condition, const ConstrExpr& constraint) {
    double coeff;
    int id;
    size_t n = 0;
    for (auto p : condition) {
	++n;
	id = p.var().id();
	coeff = p.coefficient();
    }
    if (n != 1 || !((coeff == 1 && condition.constant() == 0) || (coeff == -1 && condition.constant() == 1)))
	throw std::logic_error("condition for addConditionalConstraint must be boolean variable");
    bool complemented = coeff != 1;

    int row = CPXgetnumrows(env_, lp_);
    std::vector<int> rmatind;
    std::vector<double> rmatval;
    for (auto term : constraint.lhs()) {
	rmatind.push_back(term.var().id());
	rmatval.push_back(term.coefficient());
    }
    int nzcnt = rmatind.size();
    double rhs = -constraint.lhs().constant();
    char sense;
    switch (constraint.sense()) {
    case ConstrExpr::Sense::LE: sense = 'L'; break;
    case ConstrExpr::Sense::EQ: sense = 'E'; break;
    case ConstrExpr::Sense::GE: sense = 'G'; break;
    default: abort();
    }
    int status = CPXaddindconstr(env_, lp_, id, complemented, nzcnt, rhs, sense,
				 &rmatind[0], &rmatval[0], nullptr);
    if (status)
	error(status, "CplexModel::addConstraint");
    return makeConstraint(row);
}

void CplexModel::solve() {
    int status;
    if (CPXgetnumint(env_, lp_) || CPXgetnumbin(env_, lp_))
	status = CPXmipopt(env_, lp_);
    else
	status = CPXlpopt(env_, lp_);
    if (status)
	error(status, "CplexModel::solve");
}

Status CplexModel::status() const {
    int status = CPXgetstat(env_, lp_);
    //fprintf(stderr, "CPLEX status = %d\n", status);
    //fflush(stderr);
    switch (status) {
    case CPX_STAT_OPTIMAL:
    case CPXMIP_OPTIMAL:
    case CPXMIP_OPTIMAL_TOL:
	return Status::OPTIMAL; break;
    case CPXMIP_INFEASIBLE:
    case CPX_STAT_INFEASIBLE:
	return Status::INFEASIBLE; break;
    case CPXMIP_INForUNBD:
	return Status::INFEASIBLE_OR_UNBOUNDED; break;
    case CPXMIP_ABORT_FEAS:
    case CPXMIP_TIME_LIM_FEAS:
	return Status::ABORTED_FEASIBLE; break;
    case CPXMIP_TIME_LIM_INFEAS:
	return Status::ABORTED_NO_FEASIBLE; break;
    default:
	fprintf(stderr, "unknown CPLEX status = %d\n", status);
	abort();
	break;
    }
}

bool CplexModel::varIsIntegral(const Var& v) const {
    char xctype;
    int status = CPXgetctype(env_, lp_, &xctype, v.id(), v.id());
    if (status == CPXERR_NOT_MIP)
	return false;
    if (status)
	error(status, "CplexModel::varIsIntegral: CPXgetctype");
    switch (xctype) {
    case 'I': case 'B': return true;  break;
    case 'C':           return false; break;
    default:
	error(status, "CplexModel::varIsIntegral: unknown column type " + xctype);
	break;
    }
    return false;		// suppress warning
}

double CplexModel::varlb(const Var& v) const {
    double lb;
    int status = CPXgetlb(env_, lp_, &lb, v.id(), v.id());
    if (status)
	error(status, "CplexModel::varlbl");
    return lb;
}

double CplexModel::varub(const Var& v) const {
    double ub;
    int status = CPXgetub(env_, lp_, &ub, v.id(), v.id());
    if (status)
	error(status, "CplexModel::varub");
    return ub;
}

double CplexModel::varValue(const Var& v) const {
    double x;
    int status = CPXgetx(env_, lp_, &x, v.id(), v.id());
    if (status)
	error(status, "CplexModel::varValue");

    // CPLEX sometimes returns slightly nonintegral values even for
    // integer variables. To avoid messy logic, we properly round them
    // here.
    char xctype;
    status = CPXgetctype(env_, lp_, &xctype, v.id(), v.id());
    if (status != CPXERR_NOT_MIP) {
	switch (xctype) {
	case 'I': case 'B':
	    x = std::round(x);
	    break;
	case 'C':
	    break;
	default:
	    error(status, "CplexModel::varValue: unknown column type " + xctype);
	    break;
	}
    }

    return x;
}

double CplexModel::dualValue(const Constraint& c) const {
    double x;
    int status = CPXgetpi(env_, lp_, &x, c.id(), c.id());
    if (status)
	error(status, "CplexModel::dualValue");
    return x;
}

void CplexModel::varChangeBounds(const Var& v, double lb, double ub) const {
    const int count = 1;
    int indices[count] = { v.id() };
    double lbs[count] = { lb };
    double ubs[count] = { ub };
    int status;

    if ((status = CPXchgbds(env_, lp_, count, indices, "L", lbs)))
	error(status, "CplexModel::varChangeBounds(L)");
    if ((status = CPXchgbds(env_, lp_, count, indices, "U", ubs)))
	error(status, "CplexModel::varChangeBounds(U)");
}

void CplexModel::varSetIntegral(const Var& v, bool integral) const {
    const int indices[1] = { v.id() };
    const char xctype[1] = { integral ? CPX_INTEGER : CPX_CONTINUOUS };
    int status = CPXchgctype(env_, lp_, 1, indices, xctype);
    if (status)
	error(status, "CplexModel::varSetIntegral: CPXchgctype");
}

void CplexModel::setCoefficient(const Constraint& c, const Var& v, double x) const {
    int status = CPXchgcoef(env_, lp_, c.id(), v.id(), x);
    if (status)
	error(status, "CplexModel::setCoefficient: CPXchgcoef");
}

void CplexModel::setObjective(const Expr& obj) {
    size_t n = numVariables();
    std::vector<int> indices;
    std::vector<double> values;
    int status;
    if (n) {
	double* oldObj = new double[n];
	if ((status = CPXgetobj(env_, lp_, oldObj, 0, n - 1)))
	    error(status, "CplexModel::setObjective: CPXgetobj");
	for (size_t i = 0; i < n; i++) {
	    double newCoeff = obj.coefficient(makeVar(i));
	    if (newCoeff != oldObj[i]) {
		indices.push_back(i);
		values.push_back(newCoeff);
	    }
	}
	delete[] oldObj;
    }

    status = CPXchgobj(env_, lp_, indices.size(), &indices[0], &values[0]);
    if (status)
	error(status, "CplexModel::setObjective");
    setObjectiveConstant(obj.constant());
}
void CplexModel::setObjectiveCoeff(const Var& v, double c) {
    int indices[1] = { v.id() };
    double values[1] = { c };
    int status = CPXchgobj(env_, lp_, 1, indices, values);
    if (status)
	error(status, "CplexModel::setObjectiveCoeff");
}

void CplexModel::setObjectiveConstant(double x) {
    objConstant_ = x;
}

void CplexModel::setObjectiveSense(Sense s) {
#if CPX_VERSION_VERSION > 12 || (CPX_VERSION_VERSION == 12 && CPX_VERSION_RELEASE >= 5)
    int status = CPXchgobjsen(env_, lp_, s == Sense::MINIMIZE ? CPX_MIN : CPX_MAX);
    if (status)
	error(status, "CplexModel::setObjectiveSense");
#else
    CPXchgobjsen(env_, lp_, s == Sense::MINIMIZE ? CPX_MIN : CPX_MAX);
#endif
}

void CplexModel::setEmphasis(Emphasis e) {
    int emph = CPX_MIPEMPHASIS_BALANCED;
    switch (e) {
    case Emphasis::BALANCED:         emph = CPX_MIPEMPHASIS_BALANCED;    break;
    case Emphasis::FIND_FEASIBLE:    emph = CPX_MIPEMPHASIS_FEASIBILITY; break;
    case Emphasis::PROVE_OPTIMALITY: emph = CPX_MIPEMPHASIS_OPTIMALITY;  break;
    case Emphasis::MOVE_BOUND:       emph = CPX_MIPEMPHASIS_BESTBOUND;   break;
    }
    int status = CPXsetintparam(env_, CPXPARAM_Emphasis_MIP, emph);
    if (status)
	error(status, "CplexModel::setEmphasis");
}

void CplexModel::MIPstart(const std::map<Var, double>& vars) {
    std::vector<int> indices;
    std::vector<double> values;
    for (auto p : vars) {
        indices.push_back(p.first.id());
        values.push_back(p.second);
    }
    int ms_beg[1] = {0};
    int effortlevel[1] = {CPX_MIPSTART_REPAIR};
    CHECKSTATUS(CPXaddmipstarts, (env_, lp_, 1, indices.size(), ms_beg,
				  indices.data(), values.data(), effortlevel, NULL));
}

static void cplexError(CPXCENVptr env, int status, const std::string& msg) {
    char buffer[CPXMESSAGEBUFSIZE];
    CPXgeterrorstring(env, status, buffer);
    throw std::runtime_error(msg + std::string(": ") + std::string(buffer));
}

void CplexModel::error(int status, const std::string& msg) const {
    cplexError(env_, status, msg);
}

class CplexCallbackSolver : public CallbackSolver {
public:
    CplexCallbackSolver(CPXCENVptr env, void* cbdata, int wherefrom) :
	env_(env), cbdata_(cbdata), wherefrom_(wherefrom) {
	CHECKSTATUS(CPXgetcallbackinfo, (env, cbdata, wherefrom, CPX_CALLBACK_INFO_USER_PROBLEM, &lp_));
	std::size_t n = CPXgetnumcols(env, lp_);
	if (n == 0)
	    return;
	x_.resize(n);
	CHECKSTATUS(CPXgetcallbacknodex, (env, cbdata, wherefrom, &x_[0], 0, n - 1));
	if (wherefrom == CPX_CALLBACK_MIP_CUT_FEAS || wherefrom == CPX_CALLBACK_MIP_CUT_UNBD) {
	    // CPLEX sometimes returns slightly nonintegral values even
	    // for integer variables. To avoid messy logic in lazy
	    // constraint callbacks, we properly round them here.
	    std::vector<char> colType(n);
	    int status = CPXgetctype(env, lp_, &colType[0], 0, n - 1);
	    if (status == CPXERR_NOT_MIP)
		return;
	    if (status)
		error(status, std::string(__func__) + ": CPXgetctype");
	    for (std::size_t i = 0; i < n; ++i) {
		switch (colType[i]) {
		case 'C': break;
		case 'I': case 'B': x_[i] = intRound(x_[i]); break;
		default: throw std::runtime_error(std::string(__func__) + ": unknown column type"); break;
		}
	    }
	}
    }
    double value(const Var& v) const { return x_[v.id()]; }
    double bestBound() const {
	double bound;
	CHECKSTATUS(CPXgetcallbackinfo, (env_, cbdata_, wherefrom_,
					 CPX_CALLBACK_INFO_BEST_REMAINING, &bound));
	return bound;
    }
    double bestObjectiveValue() const {
	int haveFeasible;
	CHECKSTATUS(CPXgetcallbackinfo, (env_, cbdata_, wherefrom_,
					 CPX_CALLBACK_INFO_MIP_FEAS, &haveFeasible));
	if (!haveFeasible)
	    return CPXgetobjsen(env_, lp_) == CPX_MIN ? inf : -inf;
	double objval;
	CHECKSTATUS(CPXgetcallbackinfo, (env_, cbdata_, wherefrom_,
					 CPX_CALLBACK_INFO_BEST_INTEGER, &objval));
	return objval;
    }
    // TODO: try other parameters instead of CPX_USECUT_FORCE
    void addLazyConstraint(const ConstrExpr& constr) { addConstraint(constr, CPX_USECUT_FORCE); }
    //void addLazyConstraint(const ConstrExpr& constr) { addConstraint(constr, CPX_USECUT_PURGE); }
    void addUserCut(const ConstrExpr& constr) { addConstraint(constr, CPX_USECUT_PURGE); }

private:
    void addConstraint(const ConstrExpr& constr, int purgeable) {
	std::vector<int> indices;
	std::vector<double> values;
	for (auto term : constr.lhs()) {
	    indices.push_back(term.var().id());
	    values.push_back(term.coefficient());
	}
	double rhs = -constr.lhs().constant();
        CHECKSTATUS(CPXcutcallbackadd, (env_, cbdata_, wherefrom_, indices.size(), rhs,
					cplexConstrSense(constr.sense()),
                                        &indices[0], &values[0], purgeable));
    }

    void error(int status, const std::string& msg) const { cplexError(env_, status, msg); }
    CPXCENVptr env_;
    CPXLPptr lp_;
    void* cbdata_;
    int wherefrom_;
    std::vector<double> x_;
};

static int cplexCallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p) {
    Callback& callback = *reinterpret_cast<Callback*>(cbhandle);
    CplexCallbackSolver callbackSolver(env, cbdata, wherefrom);
    CallbackModel callbackModel(callbackSolver);
    callback.call(callbackModel);
    *useraction_p = CPX_CALLBACK_DEFAULT;
    return 0;
}

void CplexModel::setCallback(Callback& callback, CallbackType type) {
    // the callback wants to use variable indices of the original and not the presolved ILP
    CHECKSTATUS(CPXsetintparam, (env_, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF));
    if (type & USER_CUTS)
        CHECKSTATUS(CPXsetusercutcallbackfunc, (env_, cplexCallback, &callback));
    if (type & LAZY_CONSTRAINTS)
        CHECKSTATUS(CPXsetlazyconstraintcallbackfunc, (env_, cplexCallback, &callback));
#ifndef LP_SINGLE_THREAD
    // by default, CPLEX resorts to one thread in the presence of control callbacks
    // there seems to be no way to get CPLEX to determine a reasonable number of CPUs by itself
    CHECKSTATUS(CPXsetintparam, (env_, CPXPARAM_Threads, sysconf(_SC_NPROCESSORS_ONLN)));
#endif
}

void CplexModel::setTimeLimit(double seconds) {
    CHECKSTATUS(CPXsetdblparam, (env_, CPX_PARAM_TILIM, seconds));
}

} // namespace LP
