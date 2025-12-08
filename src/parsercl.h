/*************************************************************************
 ****                                                                 ****
 ****                       Include file for PARSER.CPP               ****
 ****                   Copyright 1992 by Bijarne Stroustrup          ****
 ****                    Modified 1998-2002 by Davide Bucci           ****
 ****                                                                 ****
 *************************************************************************/


/* This file is part of cRCWA.

    cRCWA is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.
    
    cRCWA is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
    details.
    
    You should have received a copy of the GNU General Public License along
    with cRCWA. If not, see <https://www.gnu.org/licenses/>. 

    Davide Bucci, 2008-2025
    Jérôme Michallon, 2012-2014
*/

#ifndef _PARSECL
#define _PARSECL


#define  PARSE_ERR_ALL            0   /* Stop on every error of any nature   */
#define  PARSE_ERR_IGNOREZERO     1   /* Ignore division by zero errors      */
#define  PARSE_ERR_IGNOREDOMAIN   2   /* Ignore zero and domain errors       */
#define  PARSE_ERR_IGNOREBADTOKEN 4   /* Ignore bad tokens in string         */
#define  PARSE_ERR_IGNORERP       8   /* Ignore ") expected" errors          */
#define  PARSE_ERR_IGNORENOTFOUND 16  /* Ignore "variable not found" error   */
#define  PARSE_ERR_IGNOREALL      32  /* Ignore all errors that may occourr  */
#define  PARSE_ERR_IGNOREOVERFLOW 128 /* Ignore all overflow errors          */
#define  PARSE_ERR_SHOWERR        64  /* Show the message of the error       */

        // Error constants
#define  ERR_ZERO            1      /* Zero divide                           */
#define  ERR_DOMAIN          2      /* Domain error sqrt(-1)                 */
#define  ERR_BADTOKEN        3      /* Unrecognized token                    */
#define  ERR_RP              4      /* ")" expected                          */
#define  ERR_NOTFOUND        5      /* Variable to parse not found           */
#define  ERR_PRIMARYEXPECTED 6      /* Expectd primary (variable or constant)*/
#define  ERR_OVERFLOW        7      /* Overflow                              */
#define  ERR_OUTOFMEMORY     8      /* Out of memory                         */
#define  ERR_LP              9      /* "(" expected                          */
#define  ERR_FLOATING        10     /* Floating point exception              */

// Questa linea è fondamentale con lo gcc
//#ifndef __TURBOC__
#define far
//#endif

#ifndef RC_INVOKED              /* Invoked by the Resource Compiler        */


#define TBLSZ 23

#define NO_DECIMAL -1


#include <stdio.h>
/*#define DEBUG(nome)   \
//    printf("Esecuzione di %s: situazione di Derive\n",nome);\
//    printf("pcExpression: %s\n", pcExpression);               \
        printf("curr_tok:     %d\n", curr_tok);                   \
        printf("no_of_errors: %d\n", no_of_errors);               \
        printf("iCurrChar:    %d\n", iCurrChar);                  \
        printf("number_value: %f\n", number_value);               \
        printf("name_string:  %s\n", name_string)
 */

struct name {
  char  * sstring;
  name  * next;
  double value;
  bool isvector;
  long int boundary;
  double *pvector;
};

class parseTok {
  int seekIndex;                // Used by GetNextVariable and ResetVariableSeek
  name * seekPos;
  struct name * table[TBLSZ];


  int searchNameOnTable(char ch);// Loads into name_string the following name.
  int extractNumber(int iCharCounter);   // Position of the number

//protected:
public:
  enum {False, True=-1};
  enum token_value {          // Operators recognized (tokens)
    NAME=3,       NUMBER,       END,
    PLUS='+',   MINUS='-',    MUL='*',    DIV='/',
    PRINT=';',  ASSIGN='=',   pLP='(',    pRP=')',
    POW='^',    LOG='l',      LOG10='L',  SQR='q',
    SIN='s',    COS='c',      TAN='t',    ABS='a',
    ARCSIN='n', ARCCOS='o',   ARCTAN='r',
    SINH='b',   COSH='d',     TANH='u',
    iINT='i',   STIRLING='g', SGN='P', COMPARE='r',
    GREATER='>',LESS='<',     NEQ='@', GEQ='`',
    LEQ='$',    AND='&',      OR='|',  NOT='!',
    SIZE='z'
    };

  int no_of_errors;             // Number of errors occourred in the elaboration
  char *pcExpression;        // Expression to parse
  int iCurrChar;                // Current char being parsed
  int curr_tok;                 // Current token found
  double number_value;
  char *name_string;            // Name of the variable being parsed
  int ErrPos;                   // Error position into string
  int ErrLen;                   // Length of the error (used in mispelled names)
  long int arrayIndex;          // Index of the array.

  int getToken(void);           // Charge a token
  name * look(const char far* p, int ins =0);
  name * insert(const char far* s) { return look(s,1);}

  int error(int);
public:
  parseTok();
  ~parseTok();
  int getErrorState(void) { return no_of_errors; }
  int getErrPos(void) { return ErrPos; }
  int getErrLen(void) { return ErrLen; }
  int setError(int ErrNo) {return error(ErrNo); }
  void resetVariableSeek(void);
  int getNextVariableName(char far* Buffer, unsigned int iBufferLen);
  int killVariable(const char far* pVarName);
  int insertVar(const char far*s, double value)
  { look(s,1)->value=value; look(s,0)->isvector=false;   return 0;}
  int insertArray(const char far*s, int size, double *values);

};

class numParser: public parseTok{
  int ErrorDetectLevel;     // Level of error detection
  double expr();
  double term();
  double sec();
  double prim();
  int error (int);

public:
  int calculate(char *pExpr, double &Result);
  inline void setErrorDetectLevel(int level)
  { ErrorDetectLevel=level; }
};


#endif
#endif

