    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

    Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include "Grid/Grid.h"

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

namespace Grid { 
  namespace QCD { 


class HmcRunner : public NerscHmcRunner {
public:

  void BuildTheAction (int argc, char **argv)

  {
    typedef WilsonImplR ImplPolicy;
    typedef MobiusFermionR FermionAction;
    typedef typename FermionAction::FermionField FermionField;

    UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

    const int Ls = 12;

    FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
    FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

    // temporarily need a gauge field
    LatticeGaugeField  U(UGrid);

    // Gauge action
    SymanzikGaugeActionR Syzaction(8.0);

    // These lines are unecessary if BC are all periodic
    std::vector<Complex> boundary = {1,1,1,-1};
    FermionAction::ImplParams Params(boundary);

    Real mass = 0.05;
    Real M5 = 1.0;
    Real b = 1.5;
    Real c = 0.5;
    FermionAction FermOp(U,*FGrid,*FrbGrid, *UGrid, *UrbGrid, mass, M5, b, c, Params);
  
    ConjugateGradient<FermionField>  CG(1.0e-8,10000);

    TwoFlavourEvenOddPseudoFermionAction<ImplPolicy> Nf2a(FermOp,CG,CG);
    TwoFlavourEvenOddPseudoFermionAction<ImplPolicy> Nf2b(FermOp,CG,CG);
    TwoFlavourEvenOddPseudoFermionAction<ImplPolicy> Nf2c(FermOp,CG,CG);
    TwoFlavourEvenOddPseudoFermionAction<ImplPolicy> Nf2d(FermOp,CG,CG);
    TwoFlavourEvenOddPseudoFermionAction<ImplPolicy> Nf2e(FermOp,CG,CG);
  

    Nf2a.is_smeared=true;
    Nf2b.is_smeared=true;
    Nf2c.is_smeared=true;
    Nf2d.is_smeared=true;
    Nf2e.is_smeared=true;

    //Collect actions
    ActionLevel<LatticeGaugeField> Level1(1);
    Level1.push_back(&Nf2a);
    Level1.push_back(&Nf2b);
    Level1.push_back(&Nf2c);
    Level1.push_back(&Nf2d);
    Level1.push_back(&Nf2e);

    ActionLevel<LatticeGaugeField> Level2(4);
    Level2.push_back(&Syzaction);

    TheAction.push_back(Level1);
    TheAction.push_back(Level2);

    Run(argc,argv);
  };

};

}}

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  int threads = GridThread::GetThreads();
  std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

  HmcRunner TheHMC;
  
  TheHMC.BuildTheAction(argc,argv);

}


