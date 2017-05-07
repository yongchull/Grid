/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/DiscLoop.hpp

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

#ifndef Hadrons_DiscLoop_hpp_
#define Hadrons_DiscLoop_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                DiscLoop                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DiscLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DiscLoopPar,
                                    std::string,    q_loop,
                                    Gamma::Algebra, gamma,
                                    std::string,    output);
};

template <typename FImpl>
class TDiscLoop: public Module<DiscLoopPar>
{
    TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TDiscLoop(const std::string name);
    // destructor
    virtual ~TDiscLoop(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(DiscLoop, TDiscLoop<FIMPL>, MContraction);

/******************************************************************************
 *                       TDiscLoop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDiscLoop<FImpl>::TDiscLoop(const std::string name)
: Module<DiscLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDiscLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q_loop};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TDiscLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiscLoop<FImpl>::setup(void)
{
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDiscLoop<FImpl>::execute(void)
{
    LOG(Message) << "Computing disconnected loop contraction '" << getName() 
                 << "' using '" << par().q_loop << "' with " << par().gamma 
                 << " insertion." << std::endl;

    CorrWriter            writer(par().output);
    PropagatorField       &q_loop = *env().template getObject<PropagatorField>(par().q_loop);
    LatticeComplex        c(env().getGrid());
    Gamma                 gamma(par().gamma);
    std::vector<TComplex> buf;
    Result                result;

    c = trace(gamma*q_loop);
    sliceSum(c, buf, Tp);

    result.gamma = par().gamma;
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf[t]);
    }

    write(writer, "disc", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_DiscLoop_hpp_
