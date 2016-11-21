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

#ifndef LP_HH_INCLUDED
#define LP_HH_INCLUDED

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace LP {

static const double inf = std::numeric_limits<double>::infinity();

enum class Solver {
    CPLEX,
    GUROBI,
};

class Model;
class Callback;

Model* newModel();
Model* newModel(Solver solver);
Model* newModel(const std::vector<Solver>& solvers);

class Expr;

class Var {
    friend Model;
    friend Expr;
public:
    Var() : m_(nullptr), v_(0) { } // want this to be able to have e.g. vector<Var>
    const Model& model() const { return *m_; }
    std::string name() const;
    int id() const { return v_; }
    bool operator==(const LP::Var& v2) const {
	if (&model() != &v2.model())
	    throw std::logic_error("cannot compare variables from different models");
	return id() == v2.id();
    }
    bool operator<(const LP::Var& v2) const {
	if (&model() != &v2.model())
	    throw std::logic_error("cannot compare variables from different models");
	return id() < v2.id();
    }

    bool isIntegral() const;
    bool isBool() const { return isIntegral() && lb() == 0 && ub() == 1; }
    double lb() const;
    double ub() const;

    double value() const;

    void changeBounds(double lb, double ub) const;
    void setIntegral(bool integral) const;

private:
    Var(const Model& m, int v) : m_(&m), v_(v) { }
    const Model* m_;
    int v_;
};

class Expr {
public:
    Expr() : m_(nullptr), constant_(0) { }
    Expr(const Var& v) : m_(&v.model()), coeff_{{v.id(), 1}}, constant_(0) { }

    bool isNull() const { return m_ == nullptr; }
    const Model& model() const { return *m_; }

    Expr operator-() const;
    Expr& operator+=(const Expr& e);
    Expr& operator-=(const Expr& e);
    Expr& operator+=(double c) { constant_ += c; return *this; }
    Expr& operator-=(double c) { constant_ -= c; return *this; }
    Expr& operator*=(double c);

    class Iter {
	friend Expr;
    public:
	Iter(const Model& m, std::map<int, double>::const_iterator it) : m_(m), it_(it) { }
	class Term {
	    friend Iter;
	public:
	    Var var() const { return var_; }
	    double coefficient() const { return coefficient_; }
	private:
	    Term(double coefficient, Var var) : coefficient_(coefficient), var_(var) { }
	    double coefficient_;
	    Var var_;
	};
	bool operator!=(const Iter& other) const { return it_ != other.it_; }
	const Iter& operator++() { ++it_; return *this; }
	Term operator*() { return Term(it_->second, Var(m_, it_->first)); }
    private:
	const Model& m_;
	std::map<int, double>::const_iterator it_;
    };
    Iter begin() const { return Iter(model(), coeff_.cbegin()); }
    Iter end() const { return Iter(model(), coeff_.cend()); }
    double coefficient(Var var) const {
	if (isNull())
	    return 0;
	if (&var.model() != m_)
	    throw std::logic_error("cannot access coefficient from different model");
	auto it = coeff_.find(var.id());
	if (it == coeff_.end())
	    return 0;
	else
	    return it->second;
    }
    size_t numCoefficients() const { return coeff_.size(); }
    double constant() const { return constant_; }
    double value() const;

    LP::Var extractVar() const {
	if (constant() || coeff_.size() != 1 || coeff_.begin()->second != 1)
	    throw std::logic_error("not variable");
	Var v(*m_, coeff_.begin()->first);
	return v;
    }

    std::pair<bool, LP::Var> extractBoolVar() const {
	if (coeff_.size() != 1)
	    throw std::logic_error("not boolean variable");
	Var v(*m_, coeff_.begin()->first);
	if (!v.isBool())
	    throw std::logic_error("not boolean variable");
	double c = coeff_.begin()->second;
	if (!((c == 1 && constant() == 0) || (c == -1 &&  constant() == 1)))
	    throw std::logic_error("not boolean variable");
	return {c == 1, v};
    }

private:
    const Model* m_;
    std::map<int, double> coeff_;
    double constant_;
};

class ConstrExpr {
public:
    enum class Sense { LE, EQ, GE };
    ConstrExpr(const Expr& e, Sense sense) : e_(e), sense_(sense) { }
    const Expr& lhs() const { return e_; }
    // note: RHS is 0, not the constant in e_.
    Sense sense() const { return sense_; }

private:
    Expr e_;
    Sense sense_;
};

class Constraint {
    friend Model;
public:
    const Model& model() const { return m_; }
    int id() const { return row_; }
    ConstrExpr constrExpr() const;
    double dualValue() const;
    void setCoefficient(const Var& v, double x) const;

private:
    Constraint(const Model& m, int row) : m_(m), row_(row) { }
    const Model& m_;
    int row_;
};

enum class Status {
    OPTIMAL,
    INFEASIBLE,
    INFEASIBLE_OR_UNBOUNDED,
    ABORTED_FEASIBLE,
    ABORTED_NO_FEASIBLE,
};

enum class Emphasis {
    BALANCED,
    FIND_FEASIBLE,
    PROVE_OPTIMALITY,
    MOVE_BOUND,
};

enum CallbackType {
    LAZY_CONSTRAINTS = 1 << 0,
    USER_CUTS        = 1 << 1,
};

class Model {
    friend Var;
public:
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    virtual ~Model() { };
    enum class Sense { MINIMIZE, MAXIMIZE };

    virtual void setVerbose(bool verbose = true) = 0;

    virtual size_t numVariables() const = 0;
    virtual size_t numConstraints() const = 0;
    virtual Expr objective() const = 0;
    virtual Sense objectiveSense() const = 0;
    virtual Var addVar(double lb, double ub, double obj, const std::string& name = "") = 0;
    Var addVar(double lb, double ub, const std::string& name = "") {
	return addVar(lb, ub, 0.0, name);
    }
    virtual Var addIntVar(double lb, double ub, double obj, const std::string& name = "") = 0;
    Var addIntVar(double lb, double ub, const std::string& name = "") {
	return addIntVar(lb, ub, 0.0, name);
    }
    Var addBoolVar(double obj, const std::string& name = "") {
	return addIntVar(0, 1, obj, name);
    }
    Var addBoolVar(const std::string& name = "") {
	return addBoolVar(0.0, name);
    }
    virtual Constraint addConstraint(const ConstrExpr& c) = 0;

    // Add "if condition then constraint".
    // condition must be a boolean variable or a negated boolean variable
    virtual Constraint addConditionalConstraint(const Expr& condition, const ConstrExpr& constraint) = 0;
    virtual ConstrExpr constrExpr(const Constraint& c) const = 0;

    virtual void setObjective(const Expr& obj) = 0;
    virtual void setObjectiveCoeff(const Var& v, double c) = 0;
    virtual void setObjectiveSense(Sense s) = 0;
    virtual void setObjectiveConstant(double x) = 0;
    void minimize(const Expr& obj) { setObjective(obj);	setObjectiveSense(Sense::MINIMIZE); }
    void maximize(const Expr& obj) { setObjective(obj);	setObjectiveSense(Sense::MAXIMIZE); }

    virtual void setEmphasis(Emphasis e) = 0;
    virtual void MIPstart(const std::map<Var, double>& vars) = 0;
    virtual void setCallback(Callback& callback, CallbackType type) = 0;
    virtual void setTimeLimit(double seconds) = 0;

    virtual void solve() = 0;
    virtual Status status() const = 0;
    virtual double objectiveValue() const = 0;
    virtual double bestBound() const = 0;

    class VarIter {
	friend Model;
    public:
	bool operator!=(const VarIter& other) const { return i_ != other.i_; }
	const VarIter& operator++() { ++i_; return *this; }
	Var operator*() { return Var(m_, i_); }
    private:
	VarIter(const Model& m, int i) : m_(m), i_(i) { }
	const Model& m_;
	int i_;
    };
    VarIter varsBegin() const { return VarIter(*this, 0); }
    VarIter varsEnd()   const { return VarIter(*this, numVariables()); }

    class ConstrIter {
	friend Model;
    public:
	bool operator!=(const ConstrIter& other) const { return i_ != other.i_; }
	const ConstrIter& operator++() { ++i_; return *this; }
	Constraint operator*() { return Constraint(m_, i_); }
    private:
	ConstrIter(const Model& m, int i) : m_(m), i_(i) { }
	const Model& m_;
	int i_;
    };
    ConstrIter constraintsBegin() const { return ConstrIter(*this, 0); }
    ConstrIter constraintsEnd()   const { return ConstrIter(*this, numConstraints()); }

    virtual bool varIsIntegral(const Var& v) const = 0;
    virtual double varlb(const Var& v) const = 0;
    virtual double varub(const Var& v) const = 0;
    virtual double varValue(const Var& v) const = 0;
    virtual double dualValue(const Constraint& c) const = 0;
    virtual void varChangeBounds(const Var& v, double lb, double ub) const = 0;
    virtual void varSetIntegral(const Var& v, bool integral) const = 0;
    virtual void setCoefficient(const Constraint& c, const Var& v, double x) const = 0;

protected:
    Model() { }
    Var makeVar(int v) const { return Var(*this, v); }
    Constraint makeConstraint(int col) { return Constraint(*this, col); }

private:
    virtual std::string varName(int v) const = 0;
};

inline std::string Var::name() const { return model().varName(id()); }
inline bool Var::isIntegral() const { return model().varIsIntegral(*this); }
inline double Var::lb() const { return model().varlb(*this); }
inline double Var::ub() const { return model().varub(*this); }
inline double Var::value() const { return model().varValue(*this); }
inline double Constraint::dualValue() const { return model().dualValue(*this); }
inline void Constraint::setCoefficient(const Var& v, double x) const { model().setCoefficient(*this, v, x); }
inline void Var::changeBounds(double lb, double ub) const { return model().varChangeBounds(*this, lb, ub); }
inline void Var::setIntegral(bool integral) const { return model().varSetIntegral(*this, integral); }

inline ConstrExpr Constraint::constrExpr() const { return model().constrExpr(*this); }


inline LP::Expr operator-(const LP::Var& v) {
    return -LP::Expr(v);
}

inline LP::Expr operator+(const LP::Expr& e1, const LP::Expr& e2) {
    LP::Expr e = e1;
    e += e2;
    return e;
}

inline LP::Expr operator-(const LP::Expr& e1, const LP::Expr& e2) {
    LP::Expr e = e1;
    e -= e2;
    return e;
}

inline LP::Expr operator+(const LP::Expr& e, double c) {
    LP::Expr r = e;
    r += c;
    return r;
}

inline LP::Expr operator+(double c, const LP::Expr& e) {
    return e + c;
}

inline LP::Expr operator-(const LP::Expr& e, double c) {
    LP::Expr r = e;
    r -= c;
    return r;
}

inline LP::Expr operator-(double c, const LP::Expr& e) {
    LP::Expr r = -e;
    r += c;
    return r;
}

inline LP::Expr operator*(double c, const LP::Expr& e) {
    LP::Expr r = e;
    r *= c;
    return r;
}

inline LP::Expr operator!(const LP::Var& v) {
    if (!v.isBool())
	throw std::logic_error("operator \"!\" only valid for boolean variables");
    return 1 - v;
}

inline LP::Expr operator!(const LP::Expr& e) {
    e.extractBoolVar();
    return 1 - e;
}

inline LP::ConstrExpr operator<=(const Expr& e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::LE);
}

inline LP::ConstrExpr operator>=(const Expr& e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::GE);
}

inline LP::ConstrExpr operator==(const Expr& e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::EQ);
}

inline LP::ConstrExpr operator<=(const Expr& e1, double e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::LE);
}

inline LP::ConstrExpr operator>=(const Expr& e1, double e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::GE);
}

inline LP::ConstrExpr operator==(const Expr& e1, double e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::EQ);
}

inline LP::ConstrExpr operator<=(double e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::LE);
}

inline LP::ConstrExpr operator>=(double e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::GE);
}

inline LP::ConstrExpr operator==(double e1, const Expr& e2) {
    return LP::ConstrExpr(e1 - e2, LP::ConstrExpr::Sense::EQ);
}

inline std::ostream& operator<<(std::ostream& out, const LP::Var& v) {
    return out << v.name();
}

class CallbackSolver {
public:
    virtual ~CallbackSolver() { }
    virtual double value(const Var& v) const = 0;
    virtual double bestBound() const = 0;
    virtual double bestObjectiveValue() const = 0;
    virtual void addLazyConstraint(const ConstrExpr& constr) = 0;
    virtual void addUserCut       (const ConstrExpr& constr) = 0;
};

class CallbackModel {
public:
    CallbackModel(CallbackSolver& callbackSolver) : callbackSolver_(callbackSolver) { }
    double value(const Var& v) const { return callbackSolver_.value(v); }
    double bestBound() const { return callbackSolver_.bestBound(); }
    double bestObjectiveValue() const { return callbackSolver_.bestObjectiveValue(); }
    void addLazyConstraint(const ConstrExpr& constr) { callbackSolver_.addLazyConstraint(constr); }
    void addUserCut       (const ConstrExpr& constr) { callbackSolver_.addUserCut(constr); }

private:
    CallbackSolver& callbackSolver_;
};

class Callback {
public:
    virtual ~Callback() { }
    virtual void call(CallbackModel& callbackModel) = 0;
};

std::ostream& operator<<(std::ostream& out, const LP::Expr& expr);
std::ostream& operator<<(std::ostream& out, const LP::ConstrExpr& constr);
std::ostream& operator<<(std::ostream& out, const LP::Model& model);

inline double intRound(double x) {
    x = std::round(x);
    if (x == 0 && std::signbit(x)) // -0 -> 0
        x = 0;
    return x;
}

}

#endif // LP_HH_INCLUDED
