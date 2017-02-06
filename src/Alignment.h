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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Utils.h"

#include <vector>
#include <string>

namespace foreseqs {

typedef struct {
  int start;
  int end;
  DataType dataType;
} PartitionDesc;

class Alignment {
public:
  Alignment(const std::string & filename);
  ~Alignment( void );

  int loadPartitionsFile(const std::string & filename);

  unsigned int getSequenceCount( void ) const;
  unsigned int getSequenceLength( int partId = -1 ) const;

  unsigned int getNumberOfPartitions( void ) { return partitionDescriptors.size(); }
  unsigned int getStartPosition( int partId = -1 ) const;
  unsigned int getEndPosition( int partId = -1 ) const;
  DataType getDataType(int partId) const;
private:
  pll_msa_t * msa;
  unsigned int sequenceCount;
  unsigned int sequenceLength;
  std::vector< PartitionDesc > partitionDescriptors;
};

#include "libpll/pll.h"
#endif
