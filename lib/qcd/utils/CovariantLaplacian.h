/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacian.h

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#ifndef COVARIANT_LAPLACIAN_H
#define COVARIANT_LAPLACIAN_H

namespace Grid {
namespace QCD {

struct LaplacianParams : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianParams, 
                                  RealD, lo, 
                                  RealD, hi, 
                                  int,   MaxIter, 
                                  RealD, tolerance, 
                                  int,   degree, 
                                  int,   precision);
  
  // constructor 
  LaplacianParams(RealD lo      = 0.0, 
                  RealD hi      = 1.0, 
                  int maxit     = 1000,
                  RealD tol     = 1.0e-8, 
                  int degree    = 10,
                  int precision = 64)
    : lo(lo),
      hi(hi),
      MaxIter(maxit),
      tolerance(tol),
      degree(degree),
      precision(precision){};
};



////////////////////////////////////////////////////////////
// Laplacian operator L on adjoint fields
//
// phi: adjoint field
// L: D_mu^dag D_mu
//
// L phi(x) = Sum_mu [ U_mu(x)phi(x+mu)U_mu(x)^dag + 
//                     U_mu(x-mu)^dag phi(x-mu)U_mu(x-mu)
//                     -2phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

template <class Impl>
class LaplacianAdjointField: public Metric<typename Impl::Field> {
  public:
  INHERIT_GIMPL_TYPES(Impl);

  private:
  OperatorFunction<typename Impl::Field> &Solver;
  LaplacianParams param;
  MultiShiftFunction PowerHalf;    
  MultiShiftFunction PowerInvHalf;    

  //typedef typename Impl::Field::vector_object vobj;
  typedef typename Impl::Field::vector_object vobj;
  typedef CartesianStencil<vobj,vobj> Stencil;

  SimpleCompressor<vobj> compressor;
  int npoint = 8;
  std::vector<int> directions    = {0,1,2,3,0,1,2,3};  // forcing 4 dimensions
  std::vector<int> displacements = {1,1,1,1, -1,-1,-1,-1};

  Stencil laplace_stencil;


 public:

  LaplacianAdjointField(GridBase* grid, OperatorFunction<GaugeField>& S, LaplacianParams& p, const RealD k = 1.0)
      : U(2*Nd, grid), Uadj(2*Nd,grid), Solver(S), param(p), kappa(k), laplace_stencil(grid, npoint, 0, directions, displacements){
        AlgRemez remez(param.lo,param.hi,param.precision);
        std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/2)"<<std::endl;
        remez.generateApprox(param.degree,1,2);
        PowerHalf.Init(remez,param.tolerance,false);
        PowerInvHalf.Init(remez,param.tolerance,true);
        
        assert(Nd==4); // forced by the stencil
        factor = kappa / (double(4 * Nd));
      };

  void Mdir(const GaugeField&, GaugeField&, int, int){ assert(0);}
  void Mdiag(const GaugeField&, GaugeField&){ assert(0);}

  void ImportGauge(const GaugeField& _U) {
    for (int mu = 0; mu < Nd; mu++) {
      U[mu]      = PeekIndex<LorentzIndex>(_U, mu);
      U[mu+4]    = Cshift(U[mu], mu, -1);// U (x-mu)
      Uadj[mu]   = adj(U[mu]);
      Uadj[mu+4] = adj(U[mu+4]);
    }
  }

  void Mref(const GaugeField& in, GaugeField& out) {
    // in is an antihermitian matrix
    // test
    //GaugeField herm = in + adj(in);
    //std::cout << "AHermiticity: " << norm2(herm) << std::endl;

    GaugeLinkField tmp(in._grid);
    GaugeLinkField tmp2(in._grid);
    GaugeLinkField sum(in._grid);

    for (int nu = 0; nu < Nd; nu++) {
      sum = zero;
      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
      GaugeLinkField out_nu(out._grid);
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(in_nu, mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * in_nu * U[mu];
        sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_nu;
      }
      out_nu = (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
      PokeIndex<LorentzIndex>(out, out_nu, nu);
    }
  }

void M(const GaugeField& in, GaugeField& out) {
    // in is an antihermitian matrix test
    // GaugeField herm = in + adj(in);
    // std::cout << "AHermiticity: " << norm2(herm) << std::endl;

    typedef typename std::remove_const<typename std::remove_reference<decltype(in._odata[0](0))>::type>::type  internal_type;
    laplace_stencil.HaloExchange(in, compressor);
    
    PARALLEL_FOR_LOOP
    for (int i = 0; i < in._grid->oSites(); i++) {
      int permute_type;
      StencilEntry *SEup, *SEdown;
      internal_type temp2, *temp;
      internal_type sum;

      for (int nu = 0; nu < Nd; nu++) {
        sum = zero;

        for (int mu = 0; mu < Nd; mu++) {

          SEup = laplace_stencil.GetEntry(permute_type, mu, i);
          if (SEup->_is_local) {
            temp = &in._odata[SEup->_offset](nu);
            if (SEup->_permute) {
              permute(temp2, *temp, permute_type);
              sum += U[mu]._odata[i]() *   temp2 * Uadj[mu]._odata[i]() - 2.0 * in._odata[i](nu);
            }
            else {
              sum += U[mu]._odata[i]() * (*temp) * Uadj[mu]._odata[i]() - 2.0 * in._odata[i](nu);
            }
          }
          else {
            sum += U[mu]._odata[i]() * laplace_stencil.CommBuf()[SEup->_offset](nu) * Uadj[mu]._odata[i]() - 2.0 * in._odata[i](nu);
          }
          
          SEdown = laplace_stencil.GetEntry(permute_type, mu+4, i);
          if (SEdown->_is_local) {
            temp = &in._odata[SEdown->_offset](nu);
            if (SEdown->_permute) {
              permute(temp2, *temp, permute_type);
              sum += Uadj[mu+4]._odata[i]() *   temp2 * U[mu+4]._odata[i]();
            }
            else {
              sum += Uadj[mu+4]._odata[i]() * (*temp) * U[mu+4]._odata[i]();
            }
          }
          else {
            sum += Uadj[mu+4]._odata[i]() * laplace_stencil.CommBuf()[SEdown->_offset](nu) * U[mu+4]._odata[i]();
          }
        } // mu loop
        out._odata[i](nu) = (1.0 - kappa) * in._odata[i](nu) - factor * sum;
      }  // nu loop
    } // parallel for , i loop
  }



  void MDeriv(const GaugeField& in, GaugeField& der) {
    // in is anti-hermitian
    RealD factor = -kappa / (double(4 * Nd));
    
    for (int mu = 0; mu < Nd; mu++){
      GaugeLinkField der_mu(der._grid);
      der_mu = zero;
      for (int nu = 0; nu < Nd; nu++){
        GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
        der_mu += U[mu] * Cshift(in_nu, mu, 1) * adj(U[mu]) * in_nu;
      }
      // the minus sign comes by using the in_nu instead of the
      // adjoint in the last multiplication
      PokeIndex<LorentzIndex>(der,  -2.0 * factor * der_mu, mu);
    } 
  }

  // separating this temporarily
  void MDeriv(const GaugeField& left, const GaugeField& right,
              GaugeField& der) {
    // in is anti-hermitian
    RealD factor = -kappa / (double(4 * Nd));

    for (int mu = 0; mu < Nd; mu++) {
      GaugeLinkField der_mu(der._grid);
      der_mu = zero;
      for (int nu = 0; nu < Nd; nu++) {
        GaugeLinkField left_nu = PeekIndex<LorentzIndex>(left, nu);
        GaugeLinkField right_nu = PeekIndex<LorentzIndex>(right, nu);
        der_mu += U[mu] * Cshift(left_nu, mu, 1) * adj(U[mu]) * right_nu;
        der_mu += U[mu] * Cshift(right_nu, mu, 1) * adj(U[mu]) * left_nu;
      }
      PokeIndex<LorentzIndex>(der, -factor * der_mu, mu);
    }
  }

  void Minv(const GaugeField& in, GaugeField& inverted){
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    Solver(HermOp, in, inverted);
  }

  void MSquareRoot(GaugeField& P){
    GaugeField Gp(P._grid);
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    ConjugateGradientMultiShift<GaugeField> msCG(param.MaxIter,PowerHalf);
    msCG(HermOp,P,Gp);
    P = Gp; 
  }

  void MInvSquareRoot(GaugeField& P){
    GaugeField Gp(P._grid);
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(*this);
    ConjugateGradientMultiShift<GaugeField> msCG(param.MaxIter,PowerInvHalf);
    msCG(HermOp,P,Gp);
    P = Gp; 
  }



 private:
  RealD kappa;
  std::vector<GaugeLinkField> U;
  std::vector<GaugeLinkField> Uadj;
  RealD factor;
};


// This is just for debugging purposes
// not meant to be used by the final users

template <class Impl>
class LaplacianAlgebraField {
 public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef SU<Nc>::LatticeAlgebraVector AVector;

  LaplacianAlgebraField(GridBase* grid, const RealD k) : 
    U(Nd, grid), kappa(k){};

  void ImportGauge(const GaugeField& _U) {
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
    }
  }

  void Mdiag(const AVector& in, AVector& out) { assert(0); }

  void Mdir(const AVector& in, AVector& out, int dir, int disp) { assert(0); }

  // Operator with algebra vector inputs and outputs
  void M(const AVector& in, AVector& out) {
    GaugeLinkField tmp(in._grid);
    GaugeLinkField tmp2(in._grid);
    GaugeLinkField sum(in._grid);
    GaugeLinkField out_mat(in._grid);
    GaugeLinkField in_mat(in._grid);

    // Reconstruct matrix
    SU<Nc>::FundamentalLieAlgebraMatrix(in, in_mat);

    sum = zero;
    for (int mu = 0; mu < Nd; mu++) {
      tmp = U[mu] * Cshift(in_mat, mu, +1) * adj(U[mu]);
      tmp2 = adj(U[mu]) * in_mat * U[mu];
      sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_mat;
    }
    out_mat = (1.0 - kappa) * in_mat - kappa / (double(4 * Nd)) * sum;
    // Project
    SU<Nc>::projectOnAlgebra(out, out_mat);
  }

 private:
  RealD kappa;
  std::vector<GaugeLinkField> U;
};


}
}

#endif
