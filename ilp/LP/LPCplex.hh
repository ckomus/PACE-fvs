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

#ifndef LPCPLEX_HH_INCLUDED
#define LPCPLEX_HH_INCLUDED

#include "LP.hh"

extern "C" {
#include <ilcplex/cplex.h>
}

namespace LP {

class CplexModel : public Model {
public:
    CplexModel(Sense sense = Sense::MINIMIZE);
    ~CplexModel();
    void setVerbose(bool verbose);
    size_t numVariables() const;
    size_t numConstraints() const;
    Expr objective() const;
    Sense objectiveSense() const;
    Var addVar(double lb, double ub, double obj, const std::string& name);
    Var addIntVar(double lb, double ub, double obj, const std::string& name);
    Constraint addConstraint(const ConstrExpr& c);
    Constraint addConditionalConstraint(const Expr& condition, const ConstrExpr& constraint);
    ConstrExpr constrExpr(const Constraint& c) const;
    void solve();
    Status status() const;
    double objectiveValue() const;
    double bestBound() const;
    bool varIsIntegral(const Var& v) const;
    double varlb(const Var& v) const;
    double varub(const Var& v) const;
    double varValue(const Var& v) const;
    double dualValue(const Constraint& c) const;
    void varChangeBounds(const Var& v, double lb, double ub) const;
    void varSetIntegral(const Var& v, bool integral) const;
    void setCoefficient(const Constraint& c, const Var& v, double x) const;
    void setObjective(const Expr& obj);
    void setObjectiveCoeff(const Var& v, double c);
    void setObjectiveSense(Sense s);
    void setObjectiveConstant(double x);
    void setEmphasis(Emphasis e);
    void MIPstart(const std::map<Var, double>& vars);
    void setCallback(Callback& callback, CallbackType type);
    void setTimeLimit(double seconds);

private:
    std::string varName(int v) const;
    Var doAddVar(double lb, double ub, double obj, const std::string& name, const char* type);

    void error(int status, const std::string& msg) const;
    CPXENVptr env_;
    CPXLPptr lp_;
    bool minimizing_;
    double objConstant_;
};

}

#endif // LPCPLEX_HH_INCLUDED
