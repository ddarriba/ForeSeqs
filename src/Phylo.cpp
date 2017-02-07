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

#include "Phylo.h"

using namespace std;

typedef struct
{
  unsigned int numberOfTaxa;
  char ** labels;
} labelStr;

static int cbSetIndices(pll_utree_t * node, void *data)
{
  labelStr * labelsData;
  char * currentLabel;
  if (pllmod_utree_is_tip(node))
  {
    labelsData = (labelStr *) data;
    currentLabel = node->label;
    for (unsigned int i=0; i<labelsData->numberOfTaxa; i++)
    {
      if (!strcmp(labelsData->labels[i], currentLabel))
      {
        /* set node, clv and scaler indices */
        node->clv_index = node->scaler_index = node->node_index = i;
        break;
      }
    }
  }
  return PLL_SUCCESS;
}

namespace foreseqs {

  Phylo::Phylo(Alignment & alignment,
               pll_utree_t * tree) :
    _alignment(alignment), _tree(tree)
  {
    labelStr labelsData;
    size_t numberOfPartitions = alignment.getNumberOfPartitions();
    _partitions.reserve(numberOfPartitions);

    size_t numberOfTaxa = alignment.getSequenceCount();

    /* update tree according to the alignment */
    labelsData.numberOfTaxa = numberOfTaxa;
    labelsData.labels = alignment.getLabels();
    pllmod_utree_traverse_apply(tree,
                                NULL,
                                NULL,
                                &cbSetIndices,
                                &labelsData);

    for (size_t part = 0; part < numberOfPartitions; part++)
    {
      pll_partition_t * partition;
      const unsigned int * pll_map = 0;

      partition = pll_partition_create (
    			numberOfTaxa,                                        /* Tip CLVs */
    			(numberOfTaxa - 2),                                  /* Inner CLVs */
    			alignment.getDataType(part) == DT_NUCLEIC ? 4 : 20,  /* States */
          alignment.getSequenceLength(part),                   /* Seq. length */
    			1,                                                   /* Models */
    			(2 * numberOfTaxa - 3),                              /* P matrices */
    			4,                                                   /* Rate categories */
    			(numberOfTaxa - 2),                                  /* Scale buffers */
    			PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_PATTERN_TIP);

      //TODO: Throw exception
    	assert (partition);

      if (partition->states = 4)
        pll_map = pll_map_nt;
      else
        pll_map = pll_map_aa;

      for (unsigned int taxonIndex = 0; taxonIndex < numberOfTaxa; taxonIndex++)
      {
        int set_states = pll_set_tip_states (partition,
                                             taxonIndex,
                                             pll_map,
                                             alignment.getSequence(taxonIndex, part));

        if (set_states == PLL_FAILURE)
        {
          pll_partition_destroy(partition);
          partition = 0;
        }
      }
    }
  }

  Phylo::~Phylo( void )
  {
    for (pll_partition_t * partition : _partitions)
    {
      pll_partition_destroy(partition);
    }
  }
}
