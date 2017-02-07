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

#define  ASCII_SIZE                128
#define  EOS                       0x00000200

#define  SYM_CR                    1 << 0
#define  SYM_LF                    1 << 1
#define  SYM_LFCR                  1 << 2
#define  SYM_DIGIT                 1 << 3
#define  SYM_CHAR                  1 << 4
#define  SYM_SPACE                 1 << 5
#define  SYM_TAB                   1 << 6
#define  SYM_EOF                   1 << 7
#define  SYM_UNKNOWN               1 << 8
#define  SYM_DOT                   1 << 9
#define  SYM_COLON                 1 << 10
#define  SYM_OPAREN                1 << 11
#define  SYM_CPAREN                1 << 12
#define  SYM_COMMA                 1 << 13
#define  SYM_SEMICOLON             1 << 14
#define  SYM_EQUAL                 1 << 15
#define  SYM_DASH                  1 << 16
#define  SYM_SLASH                 1 << 17
#define  SYM_PLUS                  1 << 18
#define  SYM_OBRACKET              1 << 19
#define  SYM_CBRACKET              1 << 20

#define  TOKEN_NUMBER              1 << 0
#define  TOKEN_STRING              1 << 1
#define  TOKEN_EOF                 1 << 2
#define  TOKEN_WHITESPACE          1 << 3
#define  TOKEN_NEWLINE             1 << 4
#define  TOKEN_UNKNOWN             1 << 5
#define  TOKEN_COLON               1 << 6
#define  TOKEN_OPAREN              1 << 7
#define  TOKEN_CPAREN              1 << 8
#define  TOKEN_FLOAT               1 << 9
#define  TOKEN_COMMA               1 << 10
#define  TOKEN_SEMICOLON           1 << 11
#define  TOKEN_EQUAL               1 << 12
#define  TOKEN_DASH                1 << 13
#define  TOKEN_SLASH               1 << 14
#define  TOKEN_OBRACKET            1 << 15
#define  TOKEN_CBRACKET            1 << 16

#define CONSUME(x)         while (token.tokenType & (x)) token = get_token (&input);
#define NEXT_TOKEN         token = get_token (&input);

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

 typedef struct
  {
    int 	        tokenType;
    const char * lexeme;
    long         len;
  } lexToken;

  static int get_next_byte (void);
  static int get_next_symbol (void);
  static lexToken get_token (int * input);

  static const char * rawtext;
  static long rawtext_size;
  static long pos = 0;
  
namespace foreseqs {

  Alignment::Alignment(const string & filename)
  {
    msa = pll_phylip_parse_msa(filename.c_str(),
                               &sequenceLength);
    sequenceCount = (unsigned int) msa->count;
  }

  Alignment::~Alignment( void )
  {
      pll_msa_destroy(msa);
  }

  int Alignment::loadPartitionsFile(const string & filename)
  {
    return PLL_SUCCESS;
  }

  char ** Alignment::getLabels( void )
  {
    return msa->label;
  }

  char * Alignment::getSequence( unsigned int taxonId, int partId )
  {
    return msa->sequence[taxonId] + partitionDescriptors[partId].start;
  }

  unsigned int Alignment::getSequenceCount( void ) const
  {
    return sequenceCount;
  }

  unsigned int Alignment::getSequenceLength( int partId ) const
  {
    if (partId >= 0)
      return partitionDescriptors[partId].end - partitionDescriptors[partId].start + 1;
    else
      return sequenceLength;
  }

  unsigned int Alignment::getStartPosition( int partId ) const
  {
    if (partId >= 0)
      return partitionDescriptors[partId].start;
    else
      return 0;
  }

  unsigned int Alignment::getEndPosition( int partId ) const
  {
    if (partId >= 0)
      return partitionDescriptors[partId].end;
    else
      return sequenceLength;
  }

  DataType Alignment::getDataType(int partId) const
  {
    return partitionDescriptors[partId].dataType;
  }





  static int get_next_byte (void)
  {
    if (pos == rawtext_size)
     {
       ++pos;
       return (EOS);
     }

    return (rawtext[pos++]);
  }

  static int get_next_symbol (void)
  {
    int ch, sym;

    ch = get_next_byte ();

    if (ch == EOS) return (SYM_EOF);
    if (ch >= ASCII_SIZE) return (SYM_UNKNOWN);

    sym = lex_table[ch];

    if (sym == SYM_LF)
     {
       if (get_next_byte() == '\n')
        {
          sym = SYM_LFCR;
        }
       else
        {
          --pos;
        }
     }

    return sym;
  }

  static lexToken get_token (int * input)
  {
    lexToken token;
    long
      start_pos,
      isFloating = 0;

    token.lexeme = rawtext + pos - 1;
    start_pos    = pos;

    switch (*input)
     {
       case SYM_SLASH:
         token.tokenType = TOKEN_SLASH;
         *input = get_next_symbol();
         break;

       case SYM_DASH:
         token.tokenType = TOKEN_DASH;
         *input = get_next_symbol();
         break;

       case SYM_EQUAL:
         token.tokenType = TOKEN_EQUAL;
         *input = get_next_symbol();
         break;

       case SYM_SEMICOLON:
         token.tokenType = TOKEN_SEMICOLON;
         *input = get_next_symbol();
         break;

       case SYM_COMMA:
         token.tokenType = TOKEN_COMMA;
         *input = get_next_symbol();
         break;

       case SYM_COLON:
         token.tokenType = TOKEN_COLON;
         *input = get_next_symbol();
         break;

       case SYM_OPAREN:
         token.tokenType = TOKEN_OPAREN;
         *input = get_next_symbol();
         break;

       case SYM_CPAREN:
         token.tokenType = TOKEN_CPAREN;
         *input = get_next_symbol();
         break;

       case SYM_OBRACKET:
         token.tokenType = TOKEN_OBRACKET;
         *input = get_next_symbol();
         break;

       case SYM_CBRACKET:
         token.tokenType = TOKEN_CBRACKET;
         *input = get_next_symbol();
         break;

       case SYM_SPACE:
       case SYM_TAB:
         do
          {
            *input = get_next_symbol();
          } while (*input == SYM_SPACE || *input == SYM_TAB);
         token.len   = pos - start_pos;
         token.tokenType = TOKEN_WHITESPACE;
         if (*input == SYM_LFCR) --token.len;
         break;

       case SYM_DIGIT:
         do
          {
            *input = get_next_symbol();
          } while (*input == SYM_DIGIT);

         if (*input == SYM_DOT)
          {
            isFloating = 1;
            do
             {
               *input = get_next_symbol ();
             } while (*input == SYM_DIGIT);
          }

         if (*input != SYM_CHAR)
          {
            token.len   = pos - start_pos;
            if (!isFloating)
              token.tokenType = TOKEN_NUMBER;
            else
              token.tokenType = TOKEN_FLOAT;
          }
         else
          {
            /* check for E notation */
            if (rawtext[pos - 1] == 'E' || rawtext[pos - 1] == 'e')
             {
               *input = get_next_symbol ();

               if (*input == SYM_PLUS || *input == SYM_DASH || *input == SYM_DIGIT)
                {
                  do
                   {
                     *input = get_next_symbol ();
                   } while (*input == SYM_DIGIT);

                  if (*input != SYM_CHAR)
                   {
                     token.len = pos - start_pos;
                     token.tokenType = TOKEN_FLOAT;
                   }
                }
               else
                {
                  token.len = pos - start_pos;
                  token.tokenType = TOKEN_STRING;
                }
             }

            if (*input == SYM_CHAR)
             {
               do {
                 *input = get_next_symbol();
               } while (*input == SYM_CHAR || *input == SYM_DIGIT || *input == SYM_DOT);
               token.len   = pos - start_pos;
               token.tokenType = TOKEN_STRING;
             }
          }

         if (*input == SYM_LFCR) --token.len;
         break;

       case SYM_CHAR:
         do
          {
            *input = get_next_symbol();
          }
         while (*input == SYM_CHAR  ||
                *input == SYM_DIGIT ||
                *input == SYM_DASH  ||
                *input == SYM_DOT);
         token.len   = pos - start_pos;
         token.tokenType = TOKEN_STRING;
         if (*input == SYM_LFCR) --token.len;
         break;

       case SYM_EOF:
         token.tokenType = TOKEN_EOF;
         break;

       case SYM_CR:
       case SYM_LF:
       case SYM_LFCR:
         do
          {
            *input = get_next_symbol();
          } while (*input == SYM_CR || *input == SYM_LFCR || *input == SYM_LF);
         token.tokenType = TOKEN_NEWLINE;
         break;
       case SYM_UNKNOWN:
       default:
         token.tokenType = TOKEN_UNKNOWN;
         break;
     }

    return (token);
  }

  static vector<PartitionDesc> * parse_partition (int * inp)
  {
      int input;
      lexToken token;
      int lines = 0;
      int partition_id = 0;
      char * tmpchar;

      input  = *inp;

      NEXT_TOKEN

      vector<PartitionDesc> * partitions = new vector<PartitionDesc>();
      while (token.tokenType != TOKEN_EOF)
      {
          ++ lines;
          PartitionDesc pi;

          CONSUME (TOKEN_WHITESPACE | TOKEN_NEWLINE)

          /* read partition type */
          if (token.tokenType != TOKEN_STRING)
          {
              cerr << "Invalid datatype in partition " << lines << endl;
              delete partitions;
              return 0;
          }

          tmpchar = (char *) calloc((size_t)token.len+10, sizeof(char));
          strncpy (tmpchar, token.lexeme, (size_t)token.len);
          tmpchar[token.len] = '\0';

          /* check first for DNA */
          if (!strcasecmp(tmpchar,"DNA") || !strcasecmp(tmpchar,"NT"))
          {
              pi.dataType   = DT_NUCLEIC;
          }
          else if (!strcasecmp(tmpchar,"PROT") || !strcasecmp(tmpchar,"AA"))
          {
              /* and  protein data */
              pi.dataType  = DT_PROTEIC;
          }
          else
          {
              cerr << "Invalid datatype in partition " << lines << ": " << tmpchar << endl;
              delete partitions;
              free (tmpchar);
              return 0;
          }
          free (tmpchar);

          NEXT_TOKEN
          CONSUME(TOKEN_WHITESPACE)

          if (token.tokenType != TOKEN_COMMA)
          {
              cerr << "Expecting ',' after datatype in partition " << lines << endl;
              delete partitions;
              return 0;
          }
          NEXT_TOKEN
                  CONSUME(TOKEN_WHITESPACE)

                  /* read partition name */
                  if (token.tokenType != TOKEN_STRING)
          {
              cerr << "Expecting partition name in partition " << lines << endl;
              delete partitions;
              return 0;
          }

          tmpchar = (char *) calloc((size_t)token.len+10, sizeof(char));
          strncpy (tmpchar, token.lexeme, (size_t)token.len);
          tmpchar[token.len] = '\0';
          pi.name = tmpchar;

          NEXT_TOKEN
                  CONSUME(TOKEN_WHITESPACE)

                  /* read equal sign */
                  if (token.tokenType != TOKEN_EQUAL)
          {
              cerr << "Expecting '=' in partition " << lines << endl;
              delete partitions;
              return 0;
          }
          NEXT_TOKEN
                  CONSUME(TOKEN_WHITESPACE)

          /* read rhs */
          if (token.tokenType != TOKEN_NUMBER)
          {
              cerr << "Invalid numerical character (partition start) in partition " << lines << endl;
              delete partitions;
              return 0;
          }
          pi.start  = (int) atoi (token.lexeme);

          NEXT_TOKEN
                  CONSUME(TOKEN_WHITESPACE)

                  if  (token.tokenType == TOKEN_DASH)
          {
              NEXT_TOKEN
                      CONSUME(TOKEN_WHITESPACE)
                      if (token.tokenType != TOKEN_NUMBER)
              {
                  cerr << "Invalid numerical character (partition end) in partition " << lines << endl;
                  delete partitions;
                  return 0;
              }
              pi.end = (int) atoi (token.lexeme);
              if (pi.end < pi.start)
              {
                  cerr << "End is smaller than Start in partition " << lines << endl;
                  delete partitions;
                  return 0;
              }
              NEXT_TOKEN
                      CONSUME(TOKEN_WHITESPACE)
          }

          CONSUME(TOKEN_WHITESPACE | TOKEN_NEWLINE)

          partitions->push_back(pi);
      }

      return (partitions);
  }
}
