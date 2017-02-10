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

#include <algorithm>

using namespace std;

typedef struct
{
  unsigned int numberOfPartitions;
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
        /* set node and clv indices */
        node->clv_index = node->node_index = i;
        break;
      }
    }
  }
  return PLL_SUCCESS;
}

static int cbLoadStoreBranches(pll_utree_t * node, void *data)
{
  unsigned int currentPartition = ((unsigned int *)data)[0];
  unsigned int nextPartition = ((unsigned int *)data)[1];
  assert(node->data);

  ((double *) node->data)[currentPartition] = node->length;
  if (nextPartition != currentPartition)
  {
    node->length = node->back->length =
      ((double *) node->data)[nextPartition];
  }

  return PLL_SUCCESS;
}

static int cbAllocate(pll_utree_t * node, void *data)
{
  unsigned int numberOfPartitions = *(unsigned int *)data;
  assert(!node->data);

  node->data = (double *) malloc(numberOfPartitions * sizeof(double));
  if (node->next)
  {
    node->next->data = node->next->next->data = node->data;
  }

  for (int i=0; i<numberOfPartitions; i++)
  {
    ((double *)node->data)[i] = node->length;
  }

  return PLL_SUCCESS;
}

namespace foreseqs {

void Phylo::updateBranches( void )
{
  unsigned int partitionData[2] = { _activePartition, _activePartition };

  pllmod_utree_traverse_apply(_tree,
                              NULL,
                              NULL,
                              &cbLoadStoreBranches,
                              partitionData);
}

void Phylo::setActivePartition( int partitionId )
{
  unsigned int partitionData[2] = { _activePartition, partitionId };

  pllmod_utree_traverse_apply(_tree,
                              NULL,
                              NULL,
                              &cbLoadStoreBranches,
                              partitionData);

  _activePartition = partitionId;
}

  Phylo::Phylo(Alignment & alignment,
               pll_utree_t * tree ,
               unsigned int numRateCategories) :
    _alignment(alignment), _tree(tree)
  {
    labelStr labelsData;
    size_t numberOfPartitions = alignment.getNumberOfPartitions();
    unsigned int j;

    _partitions.reserve(numberOfPartitions);
    _missingSequences.resize(numberOfPartitions);

    _activePartition = 0;

    size_t numberOfTaxa = alignment.getSequenceCount();

    /* update tree according to the alignment */
    labelsData.numberOfPartitions = numberOfPartitions;
    labelsData.numberOfTaxa = numberOfTaxa;
    labelsData.labels = alignment.getLabels();
    pllmod_utree_traverse_apply(tree,
                                &cbAllocate,
                                NULL,
                                &cbSetIndices,
                                &labelsData);

    cout << "Processing " << numberOfPartitions << " partitions" << endl;
    for (size_t part = 0; part < numberOfPartitions; part++)
    {
      pll_partition_t * partition;
      const unsigned int * pll_map = 0;
      unsigned int sequenceLength = alignment.getSequenceLength(part);
      partition = pll_partition_create (
    			numberOfTaxa,                                        /* Tip CLVs */
    			(numberOfTaxa - 2),                                  /* Inner CLVs */
    			alignment.getDataType(part) == DT_NUCLEIC ? 4 : 20,  /* States */
          sequenceLength,                   /* Seq. length */
    			1,                                                   /* Models */
    			(2 * numberOfTaxa - 3),                              /* P matrices */
    			numRateCategories,                                   /* Rate categories */
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
        unsigned int defSitesThreshold = 0;
        char * sequence = alignment.getSequence(taxonIndex, part);

        /* check if sequence is missing */
        unsigned char undefinedSite = '-';
        for (j = 0; j < sequenceLength; j++) {
          if (sequence[j] != undefinedSite && !defSitesThreshold--)
  					break;
        }
        if (j == sequenceLength)
        {
          cout << "Add sequence " << taxonIndex << " " << endl;
  				_missingSequences[part].push_back(taxonIndex);
        }

        int set_states = pll_set_tip_states (partition,
                                             taxonIndex,
                                             pll_map,
                                             sequence);

        if (set_states == PLL_FAILURE)
        {
          pll_partition_destroy(partition);
          partition = 0;

          //TODO Throw exception
          assert(0);
        }
      }

      _partitions.push_back(partition);
    }
    assert(_partitions.size() == numberOfPartitions);

    cout << "Go for missing branches " << endl;
    pll_utree_t ** tipNodes = (pll_utree_t **) calloc ((size_t) numberOfTaxa,
                                                       sizeof(pll_utree_t *));
    pll_utree_query_tipnodes (tree, tipNodes);

pll_utree_show_ascii(tree, PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_CLV_INDEX | PLL_UTREE_SHOW_SCALER_INDEX);
    vector< vector<pll_utree_t *> > branches;
    _missingBranches =
    Utils::findMissingBranches (tipNodes,
    														_missingSequences,
    														numberOfTaxa,
    														numberOfPartitions );

    cout << "Found " << _missingBranches[0].size() << " and " << _missingBranches[1].size() << endl;
    for (size_t part = 0; part < numberOfPartitions; part++)
    {
      for (size_t i = 0; i < _missingBranches[part].size(); i++)
      {
        cout << "  " << _missingBranches[part][i]->clv_index;
      }
      cout << endl;
    }

    //_missingBranches.inse
  }

  Phylo::~Phylo( void )
  {
    for (pll_partition_t * partition : _partitions)
    {
      if (partition)
        pll_partition_destroy(partition);
    }
  }

  double Phylo::optimizeModelParameters( int partitionId, double tolerance, double epsilon )
  {
    pll_partition_t * partition = _partitions[partitionId];
    double *empiricalFreqs, *empiricalSubstRates;
    vector<double> rateCategories;
    double logLikelihood, curLikelihood;

    unsigned int * paramsIndices;

    assert(partition);

    setActivePartition(partitionId);

    paramsIndices = (unsigned int *) calloc(partition->rate_cats, sizeof(unsigned int));
    _alpha = 1.0;

    empiricalFreqs = pllmod_msa_empirical_frequencies (partition);
    empiricalSubstRates = pllmod_msa_empirical_subst_rates (partition);

    rateCategories.resize(partition->rate_cats);
    pll_compute_gamma_cats (_alpha, partition->rate_cats, &(rateCategories[0]));

    for (int i=0; i<4; ++i)
      cout << rateCategories[i] << " ";
    cout << endl;
    for (int i=0; i<4; ++i)
      cout << empiricalFreqs[i] << " ";
    cout << endl;
for (int i=0; i<6; ++i)
  cout << empiricalSubstRates[i] << " ";
cout << endl;

    pll_set_frequencies (partition, 0, empiricalFreqs);
    pll_set_subst_params (partition, 0, empiricalSubstRates);
    pll_set_category_rates (partition, &(rateCategories[0]));

    cout << "Computing initial probability" << endl;

    logLikelihood = pllmod_utree_compute_lk(partition,
                                            _tree,
                                            paramsIndices,
                                           1,
                                           1);
    cout << "  : " << logLikelihood << endl;
    curLikelihood = logLikelihood - 10;
    while (fabs(curLikelihood - logLikelihood) > epsilon)
    {
      logLikelihood = curLikelihood;

      /* branch lengths */
      curLikelihood = -1
          * pllmod_opt_optimize_branch_lengths_iterative (partition, _tree,
                                                   paramsIndices,
                                                   1e-6, 10,
                                                   tolerance,
                                                   5,
                                                   1);

      cout << "  BL opt " << curLikelihood << endl;
      curLikelihood = -1 * pllmod_algo_opt_subst_rates (partition,
                                                   _tree,
                                                   0,
                                                   paramsIndices,
                                                   0,
                                                  0.001,
                                                  10.0,
                                                  1e7,
                                                  tolerance);
      cout << "  Rates opt " << curLikelihood << endl;

      curLikelihood = -1 * pllmod_algo_opt_alpha (partition,
                                         _tree,
                                         paramsIndices,
                                         PLLMOD_OPT_MIN_ALPHA,
                                         PLLMOD_OPT_MAX_ALPHA,
                                         &_alpha,
                                         tolerance);

      cout << "  Alpha opt " << curLikelihood << endl;

      curLikelihood = -1 * pllmod_algo_opt_frequencies (partition,
                                              _tree,
                                              0,
                                              paramsIndices,
                                              1e7,
                                              tolerance);
      cout << "  Freqs opt " << curLikelihood << endl;
    }
    logLikelihood = curLikelihood;

    updateBranches();

    cout << "Optimization done: " << logLikelihood << endl;

#ifdef PRINT_TRACE
char *newick = pll_utree_export_newick(_tree);
    cout << "PARTITION " << partitionId << ": " << newick << endl;
    free(newick);
#endif
    return logLikelihood;
  }

  bool Phylo::isBranchMissing( pll_utree_t * node, int partitionId )
  {
    if (find(_missingBranches[partitionId].begin(),
        _missingBranches[partitionId].end(), node)
          == _missingBranches[partitionId].end())
      return false;
    return true;
  }

}
