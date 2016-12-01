    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/lattice/Lattice_switch_grid.h

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

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */

#ifndef GRID_LATTICE_SWITCH_H
#define GRID_LATTICE_SWITCH_H

namespace Grid {
  template<class vobj_in, class vobj_out> inline void switch_grid(const Lattice<vobj_in>& in, Lattice<vobj_out>& out){
    GridBase* grid_in  = in._grid;
    GridBase* grid_out = out._grid;
    assert(grid_in->_gdimensions == grid_out->_gdimensions);
    assert(grid_in->_ldimensions == grid_out->_ldimensions);


PARALLEL_FOR_LOOP
    for (int l = 0; l < grid_in->lSites(); l++){
      // get the indexes
      std::vector<int> lcoor;
      grid_in->LocalIndexToLocalCoor(l,lcoor);
      int in_oSite = grid_in->oIndex(lcoor);
      int in_iSite = grid_in->iIndex(lcoor);
      int out_oSite = grid_out->oIndex(lcoor);
      int out_iSite = grid_out->iIndex(lcoor);

      // transform the internal tensors
      //out._odata[out_oSite].putlane(vobj_out::scalar_object(in._odata[in_oSite].getlane(in_iSite)),out_iSite);
      convertType(in._odata[in_oSite], out._odata[out_oSite], in_iSite, out_iSite);
    }

  }

}

#endif

