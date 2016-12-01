/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_grid_maps.cc

Copyright (C) 2015

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
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

void print_grid(GridCartesian* grid) {
  std::vector<int> oCoor;
  std::vector<int> iCoor;
  for (int lane = 0; lane < grid->Nsimd(); lane++) {
    grid->iCoorFromIindex(iCoor, lane);
    std::cout << GridLogMessage << "Lane: " << lane << "  iCoor : " << iCoor
              << std::endl;
  }
  for (int s = 0; s < grid->oSites(); s++) {
    grid->oCoorFromOindex(oCoor, s);
    std::cout << GridLogMessage << "Site: " << s << "  oCoor : " << oCoor
              << std::endl;
  }

  for (int s = 0; s < grid->lSites(); s++) {
    std::vector<int> lcoor;
    grid->LocalIndexToLocalCoor(s,lcoor);
    std::cout << GridLogMessage << "Local Site: " << s << "  LCoor : " << lcoor
              << "  oIndex : "<< grid->oIndex(lcoor)    << "  iIndex : "<< grid->iIndex(lcoor) << std::endl;
  }
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);
  GridCartesian* GridF = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexF::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGridF =
      SpaceTimeGrid::makeFourDimRedBlackGrid(GridF);

  GridCartesian* GridD = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGridD =
      SpaceTimeGrid::makeFourDimRedBlackGrid(GridD);

  GridParallelRNG pRNG_F(GridF);
  pRNG_F.SeedRandomDevice();

  std::cout << GridLogMessage << "vComplexF grid " << std::endl;
  //print_grid(GridF);

  std::cout << GridLogMessage << "vComplexD grid " << std::endl;
  //print_grid(GridD);


  LatticeGaugeField U(GridF);
  LatticeGaugeField U2(GridF);
  gaussian(pRNG_F, U);
  gaussian(pRNG_F, U2);

  auto Ulane = U._odata[0].getlane(0);



  std::cout << GridLogMessage << U._odata[0] << std::endl;

  std::cout << GridLogMessage << "Lane 0 " << Ulane << std::endl;

  std::cout << GridLogMessage << U2._odata[0] << std::endl;

  // change data and reprint
  U2._odata[0].putlane(Ulane, 0);
  std::cout << "After changing one lane" << std::endl;
  std::cout << GridLogMessage << U2._odata[0] << std::endl;

  switch_grid(U,U2);
  std::cout << "After switching fields" << std::endl;
  std::cout << GridLogMessage << U2._odata[0] << std::endl;

  // Test with different precision
  LatticeGaugeField Udouble(GridD);
  switch_grid(U,Udouble);
  std::cout << "After switching fields to double precision" << std::endl;
  std::cout << GridLogMessage << U2._odata[0] << std::endl;

  // Now lets time it to measure the overhead
  // call n times
  int ncalls = 100;
  double total_time = 0.0;
  for (int i = 0; i < ncalls; i++){
    //gaussian(pRNG_F, U);
    double t0 = usecond();
    switch_grid(U, Udouble);
    total_time += usecond() - t0;
  }

  std::cout << "Total overhead time for "<< ncalls << " calls : " << total_time << " usec" <<std::endl;
  std::cout << "Total overhead time per call "<< total_time/ncalls << " usec" << std::endl;  

}
