#include <Grid/Hadrons/Modules/MScalar/FreeProp.hpp>
#include <Grid/Hadrons/Modules/MScalar/Scalar.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MScalar;

/******************************************************************************
*                        TFreeProp implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TFreeProp::TFreeProp(const std::string name)
: Module<FreePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TFreeProp::getInput(void)
{
    std::vector<std::string> in = {par().source};
    
    return in;
}

std::vector<std::string> TFreeProp::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TFreeProp::setup(void)
{
    freeMomPropName_ = FREEMOMPROP(par().mass);
    
    if (!env().hasRegisteredObject(freeMomPropName_))
    {
        env().registerLattice<ScalarField>(freeMomPropName_);
    }
    env().registerLattice<ScalarField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TFreeProp::execute(void)
{
    ScalarField &prop   = *env().createLattice<ScalarField>(getName());
    ScalarField &source = *env().getObject<ScalarField>(par().source);
    ScalarField *freeMomProp;

    if (!env().hasCreatedObject(freeMomPropName_))
    {
        LOG(Message) << "Caching momentum space free scalar propagator"
                     << " (mass= " << par().mass << ")..." << std::endl;
        freeMomProp = env().createLattice<ScalarField>(freeMomPropName_);
        SIMPL::MomentumSpacePropagator(*freeMomProp, par().mass);
    }
    else
    {
        freeMomProp = env().getObject<ScalarField>(freeMomPropName_);
    }
    LOG(Message) << "Computing free scalar propagator..." << std::endl;
    SIMPL::FreePropagator(source, prop, *freeMomProp);
    
    if (!par().output.empty())
    {
        TextWriter            writer(par().output + "." +
                                     std::to_string(env().getTrajectory()));
        std::vector<TComplex> buf;
        std::vector<Complex>  result;
        
        sliceSum(prop, buf, Tp);
        result.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
        {
            result[t] = TensorRemove(buf[t]);
        }
        write(writer, "prop", result);
    }
}
