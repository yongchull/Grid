#ifndef Hadrons_MFermion_GaugeProp_hpp_
#define Hadrons_MFermion_GaugeProp_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugeProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TGaugeProp: public Module<GaugePropPar>
{
public:
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGaugeProp(const std::string name);
    // destructor
    virtual ~TGaugeProp(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
    SolverFn     *solver_{nullptr};
};

MODULE_REGISTER_NS(GaugeProp, TGaugeProp<FIMPL>, MFermion);

/******************************************************************************
 *                      TGaugeProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugeProp<FImpl>::TGaugeProp(const std::string name)
: Module<GaugePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    env().template registerLattice<PropagatorField>(getName());
    if (Ls_ > 1)
    {
        env().template registerLattice<PropagatorField>(getName() + "_5d", Ls_);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
    << std::endl;
    
    FermionField    source(env().getGrid(Ls_)), sol(env().getGrid(Ls_)),
    tmp(env().getGrid());
    std::string     propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    PropagatorField &prop    = *env().template createLattice<PropagatorField>(propName);
    PropagatorField &fullSrc = *env().template getObject<PropagatorField>(par().source);
    SolverFn        &solver  = *env().template getObject<SolverFn>(par().solver);
    if (Ls_ > 1)
    {
        env().template createLattice<PropagatorField>(getName());
    }
    
    LOG(Message) << "Inverting using solver '" << par().solver
    << "' on source '" << par().source << "'" << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < Nc; ++c)
    {
        LOG(Message) << "Inversion for spin= " << s << ", color= " << c
        << std::endl;
        // source conversion for 4D sources
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
                PropToFerm(source, fullSrc, s, c);
            }
            else
            {
                source = zero;
                PropToFerm(tmp, fullSrc, s, c);
                InsertSlice(tmp, source, 0, 0);
                InsertSlice(tmp, source, Ls_-1, 0);
                axpby_ssp_pplus(source, 0., source, 1., source, 0, 0);
                axpby_ssp_pminus(source, 0., source, 1., source, Ls_-1, Ls_-1);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRON_ERROR("Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm(source, fullSrc, s, c);
            }
        }
        sol = zero;
        solver(sol, source);
        FermToProp(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            PropagatorField &p4d =
            *env().template getObject<PropagatorField>(getName());
            
            axpby_ssp_pminus(sol, 0., sol, 1., sol, 0, 0);
            axpby_ssp_pplus(sol, 1., sol, 1., sol, 0, Ls_-1);
            ExtractSlice(tmp, sol, 0, 0);
            FermToProp(p4d, tmp, s, c);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
