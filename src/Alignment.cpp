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

#include "Alignment.h"

using namespace std;

static int lex_table[ASCII_SIZE] = {
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN,     SYM_TAB,      SYM_CR,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN,      SYM_LF, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*      */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/*  !"# */   SYM_SPACE, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/* $%&' */ SYM_UNKNOWN, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN,
/* ()*+ */  SYM_OPAREN,  SYM_CPAREN, SYM_UNKNOWN,      SYM_PLUS,
/* ,-./ */   SYM_COMMA,    SYM_DASH,     SYM_DOT,     SYM_SLASH,
/* 0123 */   SYM_DIGIT,   SYM_DIGIT,   SYM_DIGIT,     SYM_DIGIT,
/* 4567 */   SYM_DIGIT,   SYM_DIGIT,   SYM_DIGIT,     SYM_DIGIT,
/* 89:; */   SYM_DIGIT,   SYM_DIGIT,   SYM_COLON, SYM_SEMICOLON,
/* <=>? */ SYM_UNKNOWN,   SYM_EQUAL, SYM_UNKNOWN,      SYM_CHAR,
/* @ABC */ SYM_UNKNOWN,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* DEFG */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* HIJK */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* LMNO */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* PQRS */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* TUVW */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* XYZ[ */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,  SYM_OBRACKET,
/* \]^_ */   SYM_SLASH, SYM_CBRACKET, SYM_UNKNOWN,     SYM_CHAR,
/* `abc */ SYM_UNKNOWN,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* defg */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* hijk */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* lmno */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* pqrs */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* tuvw */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,      SYM_CHAR,
/* xyz{ */    SYM_CHAR,    SYM_CHAR,    SYM_CHAR,   SYM_UNKNOWN,
/* |}~  */    SYM_CHAR, SYM_UNKNOWN, SYM_UNKNOWN,   SYM_UNKNOWN
 };

namespace foreseqs {

  Alignment::Alignment(const & string filename)
  {
    msa = pll_phylip_parse_msa(filename.c_str(),
                               &sequenceLength);
    sequenceCount = (unsigned int) msa.count;
  }

  Alignment::~Alignment( void )
  {
      pll_msa_destroy(msa);
  }

  int Alignment::loadPartitionsFile(const & string filename)
  {
    return PLL_SUCCESS;
  }

  unsigned int Alignment::getSequenceCount( void ) const
  {
    return sequenceCount;
  }

  unsigned int Alignment::getSequenceLength( void ) const
  {
    return sequenceLength;
  }
}
