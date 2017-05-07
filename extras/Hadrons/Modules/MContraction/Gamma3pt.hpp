/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/Gamma3pt.hpp

Copyright (C) 2017

Author: Andrew Lawson    <andrew.lawson1991@gmail.com>

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
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_Gamma3pt_hpp_
#define Hadrons_Gamma3pt_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * 3pt contraction with gamma matrix insertion.
 *
 * Schematic:
 *
 *             q2           q3
 *        /----<------*------<----¬
 *       /          gamma          \
 *      /                           \
 *   i *                            * f
 *      \                          /
 *       \                        /
 *        \----------->----------/
 *                   q1
 *
 *      trace(g5*q1*adj(q2)*g5*gamma*q3)
 */

/******************************************************************************
 *                               Gamma3pt                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class Gamma3ptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Gamma3ptPar,
                                    std::string,    q1,
                                    std::string,    q2,
                                    std::string,    q3,
                                    Gamma::Algebra, gamma,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TGamma3pt: public Module<Gamma3ptPar>
{
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TGamma3pt(const std::string name);
    // destructor
    virtual ~TGamma3pt(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(Gamma3pt, ARG(TGamma3pt<FIMPL, FIMPL, FIMPL>), MContraction);

/******************************************************************************
 *                       TGamma3pt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TGamma3pt<FImpl1, FImpl2, FImpl3>::TGamma3pt(const std::string name)
: Module<Gamma3ptPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().q3};
    
    return in;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TGamma3pt<FImpl1, FImpl2, FImpl3>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TGamma3pt<FImpl1, FImpl2, FImpl3>::execute(void)
{
    LOG(Message) << "Computing 3pt contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "', '" << par().q2 << "' and '"
                 << par().q3 << "', with " << par().gamma << " insertion." 
                 << std::endl;

    CorrWriter            writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    PropagatorField3      &q3 = *env().template getObject<PropagatorField3>(par().q3);
    LatticeComplex        c(env().getGrid());
    Gamma                 g5(Gamma::Algebra::Gamma5);
    Gamma                 gamma(par().gamma);
    std::vector<TComplex> buf;
    Result                result;

    c = trace(g5*q1*adj(q2)*(g5*gamma)*q3);
    sliceSum(c, buf, Tp);

    result.gamma = par().gamma;
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }

    write(writer, "gamma3pt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_Gamma3pt_hpp_
