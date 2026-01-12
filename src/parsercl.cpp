/*************************************************************************
 ****                                                                 ****
 ****                  Main file for PARSER.CPP                       ****
 ****         Inspired by the 1992 C++ book by Bijarne Stroustrup     ****
 ****              Modified 1998-2025 by Davide Bucci                 ****
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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include <iostream>


// #define NDEBUG
#include <assert.h>


using namespace std;

#include "parsercl.h"

#ifdef __cplusplus
   typedef void (*fptr)(int);
#else
   typedef void (*fptr)();
#endif

#define ERR()  if(err_code) {error(err_code); if(err_code) return 0;}

int err_code;


/*****************************************************************************/
parseTok::parseTok(void)
{
    name_string=NULL;
    no_of_errors=0;
    curr_tok=0;
    number_value=1;
    seekPos=NULL;
    seekIndex=0;
    pcExpression = NULL;
    iCurrChar = 0;
    ErrPos =0;
    ErrLen = 0;
    arrayIndex = 0;
    for(int i=0; i<TBLSZ; ++i)
        table[i]=NULL;
}

parseTok::~parseTok()
{
    name *p=NULL;
    name *n=NULL;

    if (name_string!=NULL){
        delete[] name_string;
        name_string=NULL;
    }


    for(int i=0; i<TBLSZ; ++i) {
        if ((p=table[i])!=NULL) {
            while (p!=NULL) {
                n=p->next;
                // If an array is present, we need to scrap it.
                if (p->isvector && p->pvector!=NULL)
                    delete[] p->pvector;

                delete[] p->sstring;
                delete p;
                p=n;
            }
        }
    }

}


int parseTok::error(int iErrNo)
{
    no_of_errors = iErrNo;
    return no_of_errors;
}

int parseTok::getToken(void)
{
    char ch;
    char cNextChar;

    do {        // Eat the possibly leading spaces
        if (!(ch = pcExpression[iCurrChar++])) return curr_tok = END;
    } while (ch!='\n' && isspace(ch));

    switch (ch) {
    case ';':
    case '\n':
        return curr_tok = PRINT;

    case '!':
        if (pcExpression[iCurrChar]=='='){
            ++iCurrChar;
            return curr_tok = NEQ;
        }
        return curr_tok = NOT;
    case '*':
    case '/':
    case '+':
    case '-':
    case '(':
    case ')':
    case '^':
    case '>':
    case '<':
    case '=':
        if(ch=='=' && pcExpression[iCurrChar]=='='){
            ++iCurrChar;
            return curr_tok = COMPARE;
        } else if(ch=='>' && pcExpression[iCurrChar]=='='){
            ++iCurrChar;
            return curr_tok = GEQ;
        } else if(ch=='<' && pcExpression[iCurrChar]=='='){
            ++iCurrChar;
            return curr_tok = LEQ;
        }

        curr_tok=token_value(ch);

        return curr_tok;
    case '0':  case '1':  case '2':  case '3':  case '4':
    case '5':  case '6':  case '7':  case '8':  case '9':
    case '.':
        --iCurrChar;
        iCurrChar += extractNumber(iCurrChar);
            return curr_tok = NUMBER;

    case '&':
        if(pcExpression[iCurrChar]=='&'){
            ++iCurrChar;
            return curr_tok = AND;
        }
        return searchNameOnTable (ch);
    case '|':
        if(pcExpression[iCurrChar]=='|'){
            ++iCurrChar;
            return curr_tok = OR;
        }
        return searchNameOnTable (ch);
    case 'a':
        cNextChar = pcExpression[iCurrChar++];
        if (cNextChar=='b' && pcExpression[iCurrChar]=='s'){  // Absolute value
            ++iCurrChar;
            return curr_tok = ABS;
        }

        if (cNextChar=='r'
            && pcExpression[iCurrChar]=='c'
            && pcExpression[iCurrChar+1]=='t'
            && pcExpression[iCurrChar+2]=='a'
            && pcExpression[iCurrChar+3]=='n' ){              // Arcotangent
                iCurrChar+=4;
                return curr_tok = ARCTAN;
        }

        if (cNextChar=='r'
            && pcExpression[iCurrChar]=='c'
            && pcExpression[iCurrChar+1]=='s'
            && pcExpression[iCurrChar+2]=='i'
            && pcExpression[iCurrChar+3]=='n' ){              // Arcosine
            iCurrChar+=4;
            return curr_tok = ARCSIN;
        }

        if (cNextChar=='r'
            && pcExpression[iCurrChar]=='c'
            && pcExpression[iCurrChar+1]=='c'
            && pcExpression[iCurrChar+2]=='o'
            && pcExpression[iCurrChar+3]=='s' ){              // Arcocosine
            iCurrChar+=4;
            return curr_tok = ARCCOS;
        }

        --iCurrChar;
        return searchNameOnTable (ch);

    case 'l':
        cNextChar = pcExpression[iCurrChar++];
        if (cNextChar=='n')                       // Natural logarithm
            return curr_tok = LOG;

        if (cNextChar=='g')                       // 10-base logarithm
            return curr_tok = LOG10;

        --iCurrChar;

        return searchNameOnTable (ch);

    case 'i':
        cNextChar = pcExpression[iCurrChar++];
        if (cNextChar=='n' && pcExpression[iCurrChar]=='t'){      // Integer
            ++iCurrChar;
            return curr_tok = iINT;
        }
        --iCurrChar;

        return searchNameOnTable (ch);


    case 's':
        cNextChar = pcExpression[iCurrChar++];

        if (cNextChar=='t' && pcExpression[iCurrChar]=='r'){      // Stirling formule
            ++iCurrChar;
            return curr_tok = STIRLING;
        }

        if (cNextChar=='i' && pcExpression[iCurrChar]=='n'
             && pcExpression[iCurrChar+1]=='h'){        // Hyperbolic Sine
            ++iCurrChar;
            ++iCurrChar;
            return curr_tok = SINH;
        }

        if (cNextChar=='i' && pcExpression[iCurrChar]=='n'){      // Sine
            ++iCurrChar;
            return curr_tok = SIN;
        }

        if (cNextChar=='i' && pcExpression[iCurrChar]=='z'
             && pcExpression[iCurrChar+1]=='e'){        // Hyperbolic Sine
            ++iCurrChar;
            ++iCurrChar;
            return curr_tok = SIZE;
        }

        if (cNextChar=='q' && pcExpression[iCurrChar]=='r'){      // Square root
            ++iCurrChar;
            return curr_tok = SQR;
        }

        if (cNextChar=='g' && pcExpression[iCurrChar]=='n'){      // Sign of
            ++iCurrChar;
            return curr_tok = SGN;
        }

        --iCurrChar;

        return searchNameOnTable (ch);

    case 'c':
        cNextChar = pcExpression[iCurrChar++];

        if (cNextChar=='o' && pcExpression[iCurrChar]=='s'
            && pcExpression[iCurrChar+1]=='h'){      // Hyperbolic Cosine
            ++iCurrChar;
            ++iCurrChar;
            return curr_tok = COSH;
        }

        if (cNextChar=='o' && pcExpression[iCurrChar]=='s'){      // Cosine
            ++iCurrChar;
            return curr_tok = COS;
        }

        --iCurrChar;

        return searchNameOnTable (ch);

    case 't':
        cNextChar = pcExpression[iCurrChar++];
        if (cNextChar=='a' && pcExpression[iCurrChar]=='n'
             && pcExpression[iCurrChar+1]=='h'){      // Hyperbolic Tangent
            ++iCurrChar;
            ++iCurrChar;
            return curr_tok = TANH;
        }

        if (cNextChar=='a' && pcExpression[iCurrChar]=='n'){      // Tangent
            ++iCurrChar;
            return curr_tok = TAN;
        }

        --iCurrChar;

        return searchNameOnTable(ch);

        default:
        return searchNameOnTable(ch);
    }
}

// Loads into name_string the following name.
int parseTok::searchNameOnTable(char ch)
{
    if (name_string==NULL){
        name_string=new char[256];
        name_string[0]='\0';
        if (name_string==NULL) // Check the initialisation
            error(ERR_OUTOFMEMORY);

    }

    if (isalpha(ch))  {
        char* p = name_string;
        *p++ = ch;
        // Scanning
        while (((ch = pcExpression[iCurrChar++])!=0) && isalnum(ch)) *p++ = ch;
        --iCurrChar;
        *p = 0;       // Terminates the string
        return curr_tok=NAME;
    }
    error(ERR_BADTOKEN);

    return curr_tok=PRINT;
}

// Extract numeric constants (returns constant's length in chars)
int parseTok::extractNumber (int  iCharCounter)   // Position of the number
{
    int iStart=iCharCounter;       // Start position of the constant
    number_value=0;                 // Zero the initial value of number_value (which will contain the result)
    double mantissa=0;                // Defines the mantissa variable for Scientific-Notation constants.
    int FlagScientific=0;             // Flag: is a Scientific-Notated contant
    int FlagMinus=0;                // Flag: scientific exponent is negative
    int iPointPos=NO_DECIMAL;        // POsition of decimal separator
    char ch;               // Char read


    while (isdigit(ch=(pcExpression[iCharCounter])) ||
        (ch=='.') || (ch=='e') ||(ch=='E')) {
        switch (ch) {
            case '.':                         // Decimal separator
                ++iCharCounter;
                iPointPos =iCharCounter;
                break;
            case 'e':                   // Scientific-Notation 1E+14
            case 'E':
                FlagScientific=1;
                if (iPointPos!=NO_DECIMAL)     {
                    mantissa= number_value / pow(10, iCharCounter-iPointPos);
                    iPointPos=NO_DECIMAL;
                }else{
                    mantissa=number_value;
                }

                number_value=0;
                if (pcExpression[++iCharCounter]=='-') {
                    FlagMinus=True;
                    ++iCharCounter;
                }
                if (pcExpression[iCharCounter]=='+') {
                    ++iCharCounter;
                }
                break;
            default:
                number_value *= 10;
                number_value += int(pcExpression[iCharCounter]-'0');
                ++iCharCounter;
        }
    }

    if (FlagScientific==0) {
        if (iPointPos!=NO_DECIMAL)
            number_value /= pow(10, iCharCounter-iPointPos);
    } else {
        if (iPointPos!=NO_DECIMAL)
            mantissa /= pow(10, iCharCounter-iPointPos);

        if (FlagMinus)
            number_value= -number_value;  // Handle negative exponent
        number_value = mantissa*pow(10, number_value);
    }

    return iCharCounter-iStart;
}

// Define a new variable as an array. The contents of the array will be
// copied, so it is safe to give the pointer to an array which is
// local.
int parseTok::insertArray(const char far*s, int size, double *values)
{
    if (size<=0)
        return 0;

    name *n=look(s,1);

    // If an array is already present, we need to scrap it.
    if (n->isvector && n->pvector!=NULL)
        delete[] n->pvector;

    n->pvector=new double[size];
    n->isvector=true;
    n->boundary=size;

    for(int i=0; i<size;++i)
        n->pvector[i]=values[i];

    return 0;
}

// Search for a variable. If ins==1, eventually create it if it is not found.
name * parseTok::look(const char * p, int ins)
{
    int ii = 0;
    const char * pp=p;
    name  *n;

    assert(p!=NULL);

    // Calculate an hash value for the given string
    while (*pp) ii = ii<<1^ *pp++;
    if (ii < 0) ii = -ii;
    ii %= TBLSZ;


    for (n=table[ii]; n; n=n->next)
        if (strcmp(p,n->sstring) == 0) return n;

    if (ins == 0) {
        error(ERR_NOTFOUND);
        return NULL;
    }

    name * nn = new name;

    assert(nn!=NULL);

    nn->sstring = new char[strlen(p)+1];
    strcpy(nn->sstring, p);
    nn->value = 0;          // Assign 0 to new-created simbols
    nn->isvector=false;
    nn->pvector=NULL;
    nn->next = table[ii];
    nn->boundary=0;
    table[ii] = nn;
    return nn;
}

// Kill the given variable from variable space
int parseTok::killVariable(const char * pVarName)
{
    int ii = 0;
    const char * pp=pVarName;
    int FlagFound=False;
    name * n=NULL;
    name * Oldn=NULL;

    if(pVarName==NULL)
        return 0;

    while (*pp) ii = ii<<1^ *pp++;
    if (ii < 0) ii = -ii;
    ii %= TBLSZ;

    for (n=table[ii]; n; n=n->next) { // Seek the variable
        if (strcmp(pVarName,n->sstring) == 0) {
            FlagFound=True;
            break;
        }
                Oldn=n;     // Mantain the address of the previous variable
    }

    if (FlagFound==False) {       // Variable not found.
        return error(ERR_NOTFOUND);
    }

    if (n==table[ii])
        table[ii]=n->next;
    else
    Oldn->next=n->next;

    // Delete the variable.
    if(n->isvector && n->pvector!=NULL) {
        delete[] n->pvector;
        n->isvector=false;
    }
    delete n->sstring;
    delete n;

    return 0;
}

// Get the next variable name from variable space.
int parseTok::getNextVariableName(char * Buffer, unsigned int iBufferLen)
{
    int ret;
    int Repeat=0; // Avoid to lock the computer if the table is empty.

    while (seekPos==NULL) { // Search the next variable.
        seekPos=table[seekIndex++%TBLSZ];
        if(++Repeat>=TBLSZ+1) // Terminate the process if the table is empty,
            return 0;
    }

    strncpy(Buffer, seekPos->sstring, iBufferLen); // Copy into the buffer.
    if (strlen(seekPos->sstring)>iBufferLen)     // Calculate string lenght.
        ret=iBufferLen;
    else
        ret=strlen(seekPos->sstring);

    if (seekPos->next!=NULL)      // Seek the next variable.
        seekPos=seekPos->next;
    else {
        do {
            seekPos=table[seekIndex++%TBLSZ];
        } while (seekPos==NULL);
    }
    return ret;
}

// reset Variable seeking: the next call to GetNextVariable will return the first variable on the space.
void parseTok::resetVariableSeek(void)
{
    seekPos=NULL;
    seekIndex=0;
}

/*******************************************************************************/


int numParser::calculate(char *pcExpr, double &Result)
{
    pcExpression = pcExpr;
    iCurrChar = 0;
    no_of_errors=0;
    err_code=0;
    ErrLen=0;

    getToken();

    Result=expr();
    if (err_code) error(err_code);
    return no_of_errors;

}


double numParser::expr()
{

    double left = term();
    double right;

    static double toll=0.1;

    ERR();

    for(;;)
        switch (curr_tok) {
            case PLUS:
                getToken();
                left += term();
                break;

            case MINUS:
                getToken();
                left -= term();
                break;

            case COMPARE:
                getToken();
                right=term();
                if (left == right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case NEQ:
                getToken();
                right=term();
                if (left != right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case GREATER:
                getToken();
                right=term();
                if (left > right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case LESS:
                getToken();
                right=term();
                if (left < right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case GEQ:
                getToken();
                right=term();
                if (left >= right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case LEQ:
                getToken();
                right=term();
                if (left <= right)
                    left=1.0;
                else
                    left=0.0;
                break;

            case AND:
                getToken();
                right=term();
                if(right<0)
                    right=-right;
                if(left<0)
                    left=-left;

                if ((left>toll) && (right>toll))
                    left=1.0;
                else
                    return 0.0;

                break;

            case OR:
                getToken();
                right=term();

                if(right<0)
                    right=-right;
                if(left<0)
                    left=-left;


                if ((left>toll) || (right>toll))
                    return 1.0;
                else
                    left=0.0;
                break;

            default:
                return left;
        }
}

double numParser::term()
{
    double left;
    double d;

    left=sec();
    ERR();
    for (;;)
        switch (curr_tok) {
            case MUL:
                getToken();
                left *= sec();
                break;

            case DIV:
                getToken();
                d = sec();
                if (d==0) return error(ERR_ZERO);
                left /= d;
                break;

            default:
                return left;
        }
}

double numParser::sec()
{
    double left=prim();
    double exponent;

    ERR();
    switch (curr_tok) {

        case POW:
            getToken();
            exponent=prim();
        /*  if (exponent*left>=1000) return error(ERR_OVERFLOW);  // This is a terrible way to catch overflows! */
            if (left<0 &&  floor(exponent)!=exponent ) return error(ERR_DOMAIN);

            return pow(left, exponent);
        default:
            return left;
    }
}


double numParser::prim()
{
    double i;
    int t;
    int j;

    name * n;

    ERR();
    switch (curr_tok) {
        case NUMBER:
            getToken();
            return number_value;

        case SIZE:
            getToken();

            if(getToken()==NAME) {
                n=look(name_string, 0);

                if(n==NULL)
                    return error(ERR_NOTFOUND);

                if(!n->isvector)
                    return error(ERR_DOMAIN);
                return n->boundary;
            }

        case NAME:
            t=getToken();
            if (t == ASSIGN) { // Create the element.
                assert(name_string!=NULL);
                struct name * n = insert(name_string);
                getToken();
                n->value = expr();
                n->isvector=false;
                return n->value;
            } else if(t==pLP) { // It is an array!
                assert(name_string!=NULL);
                n=look(name_string, 1);
                // Get the index and truncate it
                double a=expr();

                arrayIndex=(long int) a;

                if (arrayIndex<0)
                    return error(ERR_DOMAIN);

                if(curr_tok == ASSIGN) {
                    getToken();
                    if(!n->isvector) {
                        n->isvector=true;
                        n->boundary=0;
                    }

                    if(arrayIndex>=n->boundary) {
                        double *np=new double[arrayIndex+1];

                        for(j=n->boundary; j<arrayIndex;++j)
                                np[j]=0;

                        if(n->isvector && n->pvector!=NULL) {
                            for(j=0; j<n->boundary;++j)
                                np[j]=n->pvector[j];

                            delete[] n->pvector;
                        }
                        n->pvector=np;
                        n->boundary=arrayIndex+1;
                    }

                    double d=expr();

                    n->pvector[arrayIndex] = d;
                } else {

                    if(!n->isvector || arrayIndex>=n->boundary)
                        return error(ERR_DOMAIN);

                    return  n->pvector[arrayIndex];
                }
            }
            if(name_string==NULL)
                return  error(ERR_NOTFOUND);
            n=look(name_string, 0);

            if (no_of_errors==ERR_NOTFOUND && ErrLen==0) {
                ErrLen=strlen(name_string);
            }
            return  n ? n->value: 0;

        case MINUS:
            getToken();

            return -term();
        case NOT:
            getToken();
            i=term();
            if(i!=0)
                return 0;
            if(i==0)
                return 1;
        case pLP:
        {   getToken();
            i = expr();
            if (curr_tok != pRP) return error (ERR_RP);
            getToken();
            return i;
        }

        case iINT:
            getToken();
            i=prim();

            if (i>=0)
                return int(i);
            else
                return int(i)-1;

        case LOG:
            getToken();
            if ((i=prim())<0)
                return error(ERR_DOMAIN);

            return log(i);

        case LOG10:
            getToken();
            if ((i=prim())<0)
                return error(ERR_DOMAIN);

            return log10(i);

        case SIN:
            getToken();
            i=prim();

            if ((i<4.6e15) && (i>-4.6e15))     // Biggest argument allowed of sin, cos, tan
                return sin(i);

            return error(ERR_OVERFLOW);

        case COS:
            getToken();
            i=prim();

            if ((i<4.6e15) && (i>-4.6e15))     // Biggest argument allowed of sin, cos, tan
                return cos(i);
            return error(ERR_OVERFLOW);

        case TAN:
            getToken();
            i=prim();

            if ((i<4.6e15) && (i>-4.6e15))     // Biggest argument allowed of sin, cos, tan
                return tan(i);
            return error(ERR_OVERFLOW);

        case SQR:
            getToken();
            if((i=prim())<0)
                return error(ERR_DOMAIN);
            else{

                return sqrt(i);
            }

        case ABS:
            getToken();

            return fabs(prim());

        case ARCSIN:
            getToken();
            i=prim();

            if (i<=1 && i>=-1)
                return asin(i);
            else
                return error(ERR_DOMAIN);

        case ARCCOS:
            getToken();
            i=prim();


            if (i<=1 && i>=-1)
                return acos(i);
            else
                return error(ERR_DOMAIN);


        case ARCTAN:
        getToken();

            return atan(prim());

        case SINH:
        getToken();

            return sinh(prim());

        case COSH:
        getToken();

            return cosh(prim());

        case TANH:
        getToken();

            return tanh(prim());

        case STIRLING:
            getToken();
            i=prim();
        if (i<0) return ERR_DOMAIN;

        return sqrt(2*M_PI*i)*pow(i,i)*pow(M_E,-i);

        case SGN:
            getToken();
        i=prim();

            if (i>0)
                return 1;
            else if (i<0)
                return -1;

            return 0;

        case END:
            return error(ERR_PRIMARYEXPECTED);

        default:
            return error(ERR_PRIMARYEXPECTED);


    }
}

int numParser::error (int iErrNo)
{

/*  ErrorDetectLevels are:
     PARSE_ERR_ALL            Stop on every error of any nature
     PARSE_ERR_IGNOREZERO     Ignore division by zero errors
     PARSE_ERR_IGNOREDOMAIN   Ignore zero and domain errors
     PARSE_ERR_INGOREBADTOKEN Ignore bad tokens in string
     PARSE_ERR_IGNORERP       Ignore ") expected" errors
     PARSE_ERR_IGNORENOTFOUND Ignore "variable not found" error
     PARSE_ERR_IGNOREALL      Ignore all errors that may occour
     PARSE_ERR_IGNOREOVERFLOW Ignore all overflow errors
     PARSE_ERR_SHOWERR        Show the message of the error
*/

    if (no_of_errors) return no_of_errors;

    switch (iErrNo) {
        case ERR_DOMAIN:    // Domain error
                        if (ErrorDetectLevel & PARSE_ERR_IGNOREDOMAIN)
                iErrNo=0;
            break;

        case ERR_ZERO:      // Zero error
            if ((ErrorDetectLevel & PARSE_ERR_IGNOREDOMAIN) ||
                (ErrorDetectLevel & PARSE_ERR_IGNOREZERO))
                iErrNo=0;
            break;

        case ERR_BADTOKEN:    // Bad token error
            if (ErrorDetectLevel & PARSE_ERR_IGNOREBADTOKEN)
                iErrNo=0;
            break;

        case ERR_RP:      // ")" expected
            if (ErrorDetectLevel & PARSE_ERR_IGNORERP)
                iErrNo=0;
            break;

        case ERR_OVERFLOW:
            if (ErrorDetectLevel & PARSE_ERR_IGNOREOVERFLOW)
                iErrNo=0;
            break;


        case ERR_NOTFOUND:    // Variable not found
            if (ErrorDetectLevel & PARSE_ERR_IGNORENOTFOUND)
                iErrNo=0;
            break;
    }

    if (ErrorDetectLevel & PARSE_ERR_IGNOREALL) iErrNo=0;

    err_code=iErrNo;

    if (iErrNo) {
    no_of_errors = iErrNo;
    ErrPos = iCurrChar;
    }
    return no_of_errors;
}

int error(int ErrNo)
{
    return err_code=ErrNo;
}
