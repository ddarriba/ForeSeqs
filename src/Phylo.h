/*
 * Model.h
 *
 *  Created on: Feb 3, 2017
 *      Author: Diego Darriba
 *      E-mail: diego.darriba@h-its.org
 *
 *  This file is part of ForeSeqs.
 *
 *  ForeSeqs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ForeSeqs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ForeSeqs.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PHYLO_H_
#define PHYLO_H_

#include "Utils.h"
#include "Alignment.h"
#include "PllDefs.h"

#include <vector>
#include <string>

namespace foreseqs {

class Phylo {
public:
  Phylo(Alignment & alignment,
        pll_utree_t * tree,
        unsigned int numRateCategories);
  ~Phylo( void );

  double optimizeModelParameters( int partitionId, double tolerance, double epsilon );
  pll_partition_t * getPartition( int partitionId ) { return _partitions[partitionId]; }
  std::vector<unsigned int> getMissingSequences( int partitionId ){ return _missingSequences[partitionId]; }
 	std::vector<pll_utree_t *> getMissingBranches( int partitionId ){ return _missingBranches[partitionId]; }

  bool isBranchMissing( pll_utree_t * node, int partitionId );
private:

  void setActivePartition( int partitionId );
  void updateBranches( void );

  Alignment & _alignment;
  std::vector< pll_partition_t * > _partitions;
  pll_utree_t * _tree;

  double _alpha;
  unsigned int _activePartition;

  std::vector<std::vector<unsigned int> > _missingSequences;
 	std::vector<std::vector<pll_utree_t *> > _missingBranches;
};

} /* namespace foreseqs */

#endif
