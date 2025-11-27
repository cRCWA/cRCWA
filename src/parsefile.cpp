/***************************************************************************
*     CLASS parsefile                                                      *
*     Davide Bucci, CROMA     2005-2013                                    *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.grenoble-inp.fr                                        *
*                                                                          *
*     Version: 1.5                                                         *
****************************************************************************
*    C++ file parsing engine                                               *
****************************************************************************/

/* Version 1.4, April 2011 - Introduce the optional "system" command */
/* Version 1.3, July 2010 - Introduce the interactive mode and include */
/* Version 1.2.2, May 2009 - Compile under GCC > 4.3 */
/* Version 1.2.1, March 2009 - Line feed styles */
/* Version 1.2, January 2009 */
/* Version 1.1, February 2006 */
/* Version 1.0, 2005 */
/* Based on ideas and code back from 1999 */

#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<cctype>
#include<string>
#include<iostream>
#include<sstream>
#include<cstring>
#include<map>
#include<stack>

using namespace std;


#include "parsercl.h"
#include "parsefile.h"


/* Pour rendre plus facile la lecture du nom de commande, on a préféré utiliser
    la même philosophie que le langage C utilise pour gérer le passage des
    paramètres entre la ligne de commande et la fonction

        int main(int argc, char *argv[])

    La variable argc contient le nombre de paramètre passés et argv[] est un
    tableau de pointeurs à chaînes de caractères, où chaque chaîne contient un
    paramètre. Le premier (d'indice 0) est le nom de la commande, comme dans
    l'exécution d'un programme lancé par la console.
    L'idée est que chaque commande est responsable de ses arguments.
*/


// Constructeur standard:
// Initialisation de la hash table
parsefile::parsefile(numParser *p)
{
    int i;
    nP=p;
    system_command_allowed=false;
    addspace=true;
    filecounter=0;
    for(i=0;i<HT_DIM;i++){
        strcpy(ht_commande[i].commande,"\0");
        ht_commande[i].fonc=NULL;
    }
    insert_command("load",this,&parsefile::c_load);
    insert_command("system",this,&parsefile::c_system);
    insert_command("label",this,&parsefile::c_label);
    insert_command("goto",this,&parsefile::c_goto);
    insert_command("for",this,&parsefile::c_for);
    insert_command("next",this,&parsefile::c_next);
    insert_command("do",this,&parsefile::c_do);
    insert_command("while",this,&parsefile::c_while);
    insert_command("if",this,&parsefile::c_if);
    insert_command("endif",this,&parsefile::c_endif);
    insert_command("else",this,&parsefile::c_else);
    insert_command("fopen",this,&parsefile::c_fopen);
    insert_command("fclose",this,&parsefile::c_fclose);
    insert_command("fscan",this,&parsefile::c_fscan);
    insert_command("fcheck",this,&parsefile::c_fcheck);
    insert_command("fskip",this,&parsefile::c_fskip);
    insert_command("fprint",this,&parsefile::c_fprint);
    insert_command("addspace",this,&parsefile::c_addspace);
}

/** LOAD: load (or include) a new file.

    Usage:
    load filename

    Parameters:
    filename complete path of the file to be included (no spaces)


*/
int parsefile::c_load(parsefile *obj, int argc,char *argv[])
{
    if(argc==2){
        if (obj->act>=MAXFILE)
            throw parsefile_commandError("load: too many files open.");

        obj->fstack[++(obj->act)] = obj->pFin = fopen(argv[1], "r");

        if (obj->pFin==NULL) {
            obj->pFin = obj->fstack[--(obj->act)];
            throw parsefile_commandError("load: can not open file.");
        }

    } else
        throw parsefile_commandError("load: invalid number of parameters.");

    return 0;
}

/** SYSTEM: execute a system command

    Usage:
    system command

    Parameters:
    command the command to be executed

    NOTE: for evident security reasons, this command is diabled by standard.
    You may explicitely activate it by calling the function
    parsefile::allow_system_command(true);
*/
int parsefile::c_system(parsefile *obj, int argc,char *argv[])
{
    if(argc>=2){
        if (!obj->system_command_allowed)
            throw parsefile_commandError("system: this command is disabled for security. Launch with the -e option.");

        // Compute the total length of the arguments
        int length=0;
        for(int i=1; i<argc; ++i) {
            length+=strlen(argv[i])+1;
        }
        char *buffer = (char *)malloc(length*sizeof(char));

        if(buffer == NULL)
            throw "system: memory error.";

        int pos = 0;
        for(int i=1; i<argc; ++i) {
            strcpy(buffer+pos, argv[i]);
            pos+=strlen(argv[i]);
            buffer[pos++]=' ';
        }

        buffer[--pos]='\0';
        int r = system(buffer);
        cout << "System command "<<buffer<<
            " executed and returned "<<r<<"."<<endl;

        free(buffer);

    } else
        throw parsefile_commandError("system: invalid number of parameters.");

    return 0;
}

/* Label */

int  parsefile::c_label(parsefile *obj, int argc,char *argv[])
{
    if(argc==2){
        long int l=obj->getBookmark();
        obj->labels[argv[1]] = l;
    } else
        throw parsefile_commandError("label: invalid number of parameters.");

    return 0;
}

/* Goto */

int  parsefile::c_goto(parsefile *obj, int argc,char *argv[])
{
    long int p=0;
    if(argc==2){
        p=obj->labels[argv[1]];
        if(p==0)
            throw parsefile_commandError("goto: label not found.");
        obj->setBookmark(p);
    } else
        throw parsefile_commandError("goto: invalid number of parameters.");

    return 0;
}

/* if

    Usage:
        if condition
            {block of code executed if condition!=0...}
        else
            {block of code executed if condition==0...}
        endif

*/
int  parsefile::c_if(parsefile *obj, int argc,char *argv[])
{
    static const double toll=0.1;   // Tolerance to consider true a condition
    if(argc==2){
        double condition=0;

        sscanf(argv[1],"%20lf", &condition);

        // Create a label, unique.
        string labelif="iflbl";
        char buffer[256];
        sprintf(buffer, "%li", obj->currentPosition);

        labelif+= buffer;

        obj->ifstack.push(labelif);
        //cout<<"if label: "<<labelif<<" p="<<obj->currentPosition<<endl;

        if(condition<0.0) condition=-condition;
        //cout<<"condition: "<<condition<<endl;

        if(condition<toll && !obj->noForCycles) {
            // jump only if the condition is "false", i.e. equal to zero.
            long int p=obj->labels[labelif];
            if(p==0) {
                // This can happen when else is not used. In this case, search
                // directly for the "endif" label
                labelif+="_e";
                p=obj->labels[labelif];

                if (p==0)
                    throw parsefile_commandError("if: can not find endif.");
            }
            obj->setBookmark(p);
        }
    } else {
        for(int i=0;i<argc;++i)
            cerr <<"parameter "<<i<<"  ["<< argv[i]<<"]"<<endl;
        throw parsefile_commandError("if: invalid number of parameters.");
    }

    return 0;
}

/* else

    Usage:
        if condition
            {block of code executed if condition!=0...}
        else
            {block of code executed if condition==0...}
        endif

*/
int  parsefile::c_else(parsefile *obj, int argc,char *argv[])
{
    if(argc==1){
        if(obj->ifstack.empty()){
            throw parsefile_commandError("else without if.");
        }
        string labelif=obj->ifstack.top();
        obj->labels[labelif] = obj->getBookmark();
        //cout<<" else saved position "<<obj->labels[labelif]<<" in labelif="<<
        //    labelif<<endl;
        if(!obj->noForCycles) {
            // Search for the next endif.
            labelif+="_e";
            long int p=obj->labels[labelif];
            // cout << "else label for endif:"<<labelif<<" p="<<p<<endl;
            if(p==0)
                throw parsefile_commandError("else: can not find endif.");
            else {
                obj->ifstack.pop();
                obj->setBookmark(p);
            }
        }

        //cout << "position: " << obj->labels[labelif]<<endl;

    } else
        throw parsefile_commandError("else: invalid number of parameters.");

    return 0;
}

/* endif

    Usage:
        if condition
            {block of code executed if condition!=0...}
        else
            {block of code executed if condition==0...}
        endif

*/
int  parsefile::c_endif(parsefile *obj, int argc,char *argv[])
{
    // NOTE: memory leak in the ifstack???

    if(argc==1){
        if(obj->ifstack.empty()){
            throw parsefile_commandError("endif without if.");
        }
        string labelif=obj->ifstack.top();
        obj->ifstack.pop();


        // Labels for endif have a _e suffix to differentiate them with those
        // used by else.
        labelif +="_e";

        obj->labels[labelif] = obj->getBookmark();
        //cout<<"endif saved position "<<obj->labels[labelif]<<" in labelif="<<
        //    labelif<<endl;

    } else
        throw parsefile_commandError("endif: invalid number of parameters.");

    return 0;
}

/* For

    Usage: for :variable start stop increment

    for(variable=start; variable<stop; variable+=increment)
*/
int  parsefile::c_for(parsefile *obj, int argc,char *argv[])
{
    if(argc==5){
        char buffer[256];
        // Get the name of the variable skipping the first character
        // (usually a ':'), even if it can be anything
        string variable=argv[1]+1;
        snprintf(buffer, 256, "%li", obj->currentPosition);

        string labelfor=variable+"+"+buffer;

        obj->forstack.push(labelfor);

        double start;
        double stop;
        double increment;

        sscanf(argv[2],"%20lf", &start);
        sscanf(argv[3],"%20lf", &stop);
        sscanf(argv[4],"%20lf", &increment);
        obj->labels[labelfor] = obj->getBookmark();

        snprintf(buffer, 256, "%s=%f", argv[1]+1, start);

        obj->nP->calculate(buffer, start);

        obj->for_end[variable] = stop;
        obj->for_increment[variable] = increment;

        if(start>=stop) {
            long int p=obj->for_next[labelfor];
            if(!obj->noForCycles) {
                obj->setBookmark(p);
                obj->forstack.pop();
            }
        }
    } else
        throw parsefile_commandError("for: invalid number of parameters.");

    return 0;
}

/* Next */
int  parsefile::c_next(parsefile *obj, int argc,char *argv[])
{
    if(argc==1 || argc==2){
        string variable;
        string labelfor=obj->forstack.top();

        int pos=labelfor.find_first_of("+",0);
        string varcheck=labelfor.substr(0,pos);

        // If there are two arguments, the second one must match the
        // variable
        if(argc==2) {
            variable = (argv[1]+1);

            if(varcheck.compare(variable)!=0) {
                string err="for: incorrect variable. It should be "+varcheck+
                    " instead of "+variable;
                throw parsefile_commandError(err.c_str());
            }
        } else {
            variable=varcheck;
        }

        double res=0;

        char buffer[256];
        snprintf(buffer, 256, "%s=%s+%f", variable.c_str(), variable.c_str(),
            obj->for_increment[variable]);

        obj->nP->calculate(buffer, res);
        obj->for_next[labelfor] = obj->getBookmark();

        if(res<obj->for_end[variable]) {
            long int p=obj->labels[labelfor];
            if(!obj->noForCycles) {
                obj->setBookmark(p);
            } else {
                obj->forstack.pop();
            }
        } else {
            obj->forstack.pop();
        }
    } else
        throw parsefile_commandError("next: invalid number of parameters.");

    return 0;
}


/** Do: start a loop (closed by while)

    Usage:
    do

*/
int  parsefile::c_do(parsefile *obj, int argc,char *argv[])
{
    if(argc==1){
        char buffer[256];

        string variable="do_var";
        snprintf(buffer, 256, "%li", obj->currentPosition);

        string labelfor=variable+"+"+buffer;
        obj->forstack.push(labelfor);
        obj->labels[labelfor] = obj->getBookmark();

    } else
        throw parsefile_commandError("do: invalid number of parameters.");

    return 0;
}

/** while: close the loop started by do.

    usage: while condition

    example

    do

    while condition
 */
int  parsefile::c_while(parsefile *obj, int argc,char *argv[])
{
    if(argc==2){
        string variable;
        string labelfor=obj->forstack.top();

        int pos=labelfor.find_first_of("+",0);
        string varcheck=labelfor.substr(0,pos);

        double res=0;

        double condition;
        sscanf(argv[1],"%20lf", &condition);

        if(condition!=0.0) {
            long int p=obj->labels[labelfor];
            if(!obj->noForCycles) {
                obj->setBookmark(p);
            } else {
                obj->forstack.pop();
            }
        } else {
            obj->forstack.pop();
        }
    } else
        throw parsefile_commandError("while: invalid number of parameters.");

    return 0;
}


/** ADDSPACE add spaces and newlines when printing

    Usage:
    addspace activate

    Parameters:
        activate select the mode, i (default) means that newlines and spaces
            are automatically added to printing commands (print and fprint).
            o means that newlines and spaces are not automatically added.
*/
int parsefile::c_addspace(parsefile *obj, int argc,char *argv[])
{
    if(argc==2){
        if(strcmp(argv[1],"i")==0) {
            obj->addspace=true;
        } else if(strcmp(argv[1],"o")==0) {
            obj->addspace=false;
        } else
            throw parsefile_commandError("addspace: unrecognized mode option. "
             "It should be {i|o}.");

    } else {
        throw parsefile_commandError("addspace: invalid number of parameters.");
    }
    return 0;
}

/** PRINT print a line at the output

    Usage:
    print wathever you want
*/
int parsefile::c_print(parsefile *obj, int argc,char *argv[])
{
    // The count starts from 1 since we do not want to see "print" at the
    // beginning of the line
    for (int i=1; i<argc; ++i) {
        cout << argv[i];
        if(obj->addspace)
            cout<< " ";
    }
    if(obj->addspace) cout << "\n";
    return 0;
}

/** FOPEN open a file

    Usage:
    fopen filename mode

    Return value: ans will contain a handle to the opened file

    Parameters:
        filename: the name of the file to be opened
        mode: the open mode. It must be "w" for writing in a file or "r" for
            reading.

    Example:
        Write a message on "test.txt" file

        fopen "test.txt" "w"
        let file=ans

        fprint file "This is a test."
        fclose file

*/
int parsefile::c_fopen(parsefile *obj, int argc,char *argv[])
{
    FILE *file;

    if(argc==3){
        double handle=obj->filecounter++;
        file=fopen(argv[1], argv[2]);

        if(file==NULL){
            throw parsefile_commandError("fopen: can not open file.");
        }
        obj->openedfiles[handle]=file;
        obj->insertVar("ans",handle);

    } else {
        throw parsefile_commandError("fopen: invalid number of parameters.");
    }
    return 0;
}

/** FCLOSE close a file

    Usage:
    fclose handle
*/
int parsefile::c_fclose(parsefile *obj, int argc,char *argv[])
{
    FILE *file;

    if(argc==2){
        double handle;
        if(sscanf(argv[1], "%20lf", &handle)!=1){
            throw parsefile_commandError("fclose: can not read file handle.");
        }

        file=obj->openedfiles[handle];

        fclose(file);

        obj->openedfiles.erase(handle);
    } else {
        throw parsefile_commandError("fclose: invalid number of parameters.");
    }
    return 0;
}

/** FSKIP skip a line on the file

    Usage:
    fskip handle
*/
int parsefile::c_fskip(parsefile *obj, int argc,char *argv[])
{
    FILE *file;

    if(argc==2){
        double handle;
        if(sscanf(argv[1], "%20lf", &handle)!=1){
            throw parsefile_commandError("fskip: can not read file handle.");
        }

        file=obj->openedfiles[handle];

        // skip a line
        fscanf(file, "%*[^\n]");

    } else {
        throw parsefile_commandError("fskip: invalid number of parameters.");
    }
    return 0;
}

/** FPRINT write on a file

    Usage:
    fprint handle "what to print"
*/
int parsefile::c_fprint(parsefile *obj, int argc,char *argv[])
{
    FILE *file;
    char space[]=" ";

    if(argc>2){
        double handle;
        if(sscanf(argv[1], "%20lf", &handle)!=1){
            throw parsefile_commandError("fprint: can not read file handle.");
        }

        file=obj->openedfiles[handle];

        if(!obj->addspace)
            space[0]='\0';

        for (int i=2; i<argc; ++i) {
            fprintf(file,"%s%s",argv[i], space);
        }
        if(obj->addspace) fprintf(file,"\n");

    } else {
        throw parsefile_commandError("fprint: invalid number of parameters.");
    }
    return 0;
}

/** FSCAN read from a file a numeric value and put it in the "ans" variable.

    Usage:
    fscan handle

    Return:
    ans contain the read value
*/
int parsefile::c_fscan(parsefile *obj, int argc,char *argv[])
{
    FILE *file;
    if(argc==2){
        double handle;
        if(sscanf(argv[1], "%20lf", &handle)!=1){
            throw parsefile_commandError("fscan: can not read file handle.");
        }

        file=obj->openedfiles[handle];

        double val=0;

        if (fscanf(file,"%lf",&val)!=1) {
            obj->lastRead=false;
        } else {
            obj->nP->insertVar("ans",val);
            obj->lastRead=true;
        }
    } else {
        throw parsefile_commandError("fscan: invalid number of parameters.");
    }
    return 0;
}

/** FCHECK check if the last read in a file has been successful

    Usage:
    fcheck

    Return:
    ans nonzero if the read is successful or 0 otherwise.
*/
int parsefile::c_fcheck(parsefile *obj, int argc,char *argv[])
{
    if(argc==1){
        double val=0;
        if(obj->lastRead)
            val=1.0;
        else
            val=0.0;

        obj->nP->insertVar("ans",val);
    } else {
        throw parsefile_commandError("fcheck: invalid number of parameters.");
    }
    return 0;
}

/** Updates the current position in the file.
*/
void parsefile::setBookmark(long int b)
{
    currentPosition=b;
    //cout << "setBookmark => "<<currentPosition<<endl;
    if(pFin!=NULL) {
        fseek(pFin, b,0);
    } else {
        throw parsefile_commandError("parsefile::setBookmark");
    }
}

/** Gets the current position in the file.
*/
long int parsefile::getBookmark(void)
{
    //cout << "getBookmark => "<<currentPosition<<endl;
    return currentPosition;
}


/**
    Read and execute the commands contained in the specified file.

    @param fileName the name of the file to be read.
    @param interactive if true, activate the interactive mode.
*/
void parsefile::read(string fileName, bool interactive)
{
    cout<<"Reading file: "<<fileName<<"\n";
    cout.flush();
    if (interactive)
        pFin = stdin;
    else {
        const char *name=fileName.c_str();
        pFin=fopen(name, "r");
    }
        // Can not open the given file.
    if (pFin==NULL){
        string errstr("File ");
        errstr+=fileName;
        errstr+=" not found";

        throw parsefile_commandError(errstr);
    }
    readFile(pFin,interactive);
}

/**
    Read and execute the commands contained in the specified file.

    @param file the file to be read.
    @param interactive if true, activate the interactive mode.
*/
void parsefile::readFile(FILE* pFin_r, bool interactive)
{
    pFin=pFin_r;
    nP->setErrorDetectLevel(0);

    char buffer[BUFDIM];
    char com[BUFDIM];
    int index,i;
    int ligne=0;
    int c_argc;
    char *c_argv[BUFDIM];
    bool to_free[BUFDIM];
    bool dont_evaluate[BUFDIM];

    // Reset some arrays
    for(i=0;i<BUFDIM;++i) {
        to_free[i] = false;
        dont_evaluate[i]=false;
    }

    bool inQuotes = false;
    fstack[0] = pFin;
    act = 0;

    if (pFin==stdin)
        cout << ligne << " > ";

    bool stopEx=false;
    // The following code is quite intricate, because two passes are to be
    // done for scripts, in order to find all the possible labels defined.
    // This is not done if the execution is stopped or if the software is
    // used in the interactive mode.
    for (int pass=0; !stopEx && pass<2; ++pass) {
        if (pFin!=stdin)
            fseek(pFin,0,0);

        bool cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
        currentPosition=ftell(pFin);

        bool doNotExecute=false;
        noForCycles=false;

        while(cont && !stopEx){
            ++ligne;
            doNotExecute=false;
            noForCycles=false;

            inQuotes = false;
            if (sscanf(buffer,"%128s",com)==1) {
                if(com[0]=='#')
                    com[1]='\0';

                index=ht_find(com);
                //cout <<"pass: "<<pass<< "; command: "<<com<<endl;
                if(index>=0){
                    // In the first pass, we only interpret and execute
                    // the following commands
                    if(pass==0&&pFin!=stdin&&
                        (strcmp(com,"label")==0 || strcmp(com,"for")==0 ||
                        strcmp(com,"next")==0 || strcmp(com,"if")==0 ||
                        strcmp(com,"endif")==0 || strcmp(com,"else")==0 ||
                        strcmp(com,"load")==0)) {
                        // Just introduce the labels in the for cycles, but do
                        // not execute them.
                        noForCycles=true;
                    } else if (pass==0&&pFin!=stdin) {
                        doNotExecute=true;
                    }

                    // Valid command: we need to separate the arguments

                    i=0;

                    // On utilise la variable suivante comme un compteur
                    // des paramètres, pour la valeur 0 on a le commande

                    c_argc=0;

                    // Every element of the array contains a parameter, but the
                    // first, which is the command itself.

                    while ((buffer[i]==' '|| buffer[i]=='\t') 
                        && i<BUFDIM-1 && !inQuotes)
                        ++i;

                    // Le prèmier élément de c_argv pointe à la commande
                    c_argv[0]=&buffer[i];

                    bool inComment=false;

                    while(buffer[i]!='\0'&&(buffer[i]!='\n' &&
                        buffer[i]!='\r')){

                        if(buffer[i]=='#') { // A comment.
                            inComment=true;
                        }
                        if(inComment) {
                            ++i;
                            continue;
                        }
                        if(buffer[i]=='"') {
                            inQuotes = !inQuotes;

                            // If a parameter does include double quotes, it
                            // will not be evaluated.

                            dont_evaluate[c_argc]=true;

                            // Make sort that the parameter passed does not
                            // include the quotes

                            if (inQuotes && i<BUFDIM-1)
                                c_argv[c_argc]=&buffer[i+1];
                            else if(!inQuotes)
                                buffer[i]='\0';
                        }
                        // on sait que chaque paramétre est separé de l'autre
                        // par un espace donc on peut chercher sur
                        // l'array "buffer" chaque paramétre

                        if((buffer[i]==' ' || buffer[i]=='\t') && !inQuotes){
                            // We substitute space with the '\0' terminator
                            // This is useful since we can make c_argv point at
                            // the same memory without having to move things
                            // around

                            buffer[i]='\0';

                            // on "mange" les espaces suivantes

                            while ((buffer[i+1]==' ' || buffer[i+1]=='\t') 
                                && i<BUFDIM-1)
                                ++i;

                            if(buffer[i+1]=='#') { // A comment.
                                break;
                            }


                            if (buffer[i+1]=='\n' || buffer[i+1]=='\0')
                                break;


                            if (++c_argc>=BUFDIM) {
                                throw parsefile_commandError("Command line too"
                                    " long.");
                            };
                            // quand on trouve l'espace, on donne l'adresse
                            // du prochaine paramètre à c_argv:

                            c_argv[c_argc]=&buffer[i+1];
                        }
                        ++i;
                        if (i>=BUFDIM){
                            throw parsefile_commandError("Command line too"
                                " long.");
                        }
                    }
                    buffer[i]='\0';

                    double res;

                    // Process all arguments (calculate, if possible)
                    for(int k=0; k<c_argc+1; ++k) {
                        int err=1;

                        if (!dont_evaluate[k]) {
                            err=nP->calculate(c_argv[k], res);
                        }
                        dont_evaluate[k]=false;
                        to_free[k] = false;

                        if (err==0) {
                            to_free[k] = true;
                            c_argv[k] = new char[BUFDIM];

                            snprintf(c_argv[k], BUFDIM, "%g", res);
                        }
                    }

                    if(!doNotExecute) {
                        try {
                            if(!interactive) {
                                if(ht_commande[index].fonc(
                                    ht_commande[index].obj,
                                    c_argc+1,c_argv)){
                                    ostringstream errstr;

                                    errstr<<"Execution error in command: ";
                                    errstr<<ht_commande[index].commande;
                                    errstr<<" at line ";
                                    errstr<<ligne;
                                    // error while executing command

                                    throw parsefile_commandError(errstr.str());
                                }
                            } else {
                                try {
                                    if(ht_commande[index].fonc(
                                        ht_commande[index].
                                        obj, c_argc+1,c_argv)){
                                        cerr<<"Execution error in command: ";
                                        cerr<<ht_commande[index].commande;
                                        cerr<<endl;
                                    }
                                } catch(parsefile_commandError P) {
                                    cerr<<P.getMess()<<"\n";
                                    if(pFin!=stdin) {
                                        cont=false;
                                        pass=1;
                                        goto test_end_file;
                                    }
                                }
                            }
                        } catch(parsefile_stop Q) {
                            cout<< "Stopping the program execution"<<endl;
                            stopEx = true;
                        }
                    }
                    for(int k=0; k<c_argc+1; ++k) {
                        if (to_free[k]) {
                            delete[] c_argv[k];
                            to_free[k]=false;
                        }
                    }
                } else {
                    cerr<<"Could not recognize the command ";
                    cerr<<com;
                    cerr<<" at line ";
                    cerr<<ligne<<"\n";

                    if(pFin!=stdin)
                        stopEx=true;
                }
            }

            if (!stopEx && (pFin==stdin))
                cout << ligne << " > ";

            // Do not worry if you do not understand.
            if (!stopEx)
                cont = (fgets(buffer, BUFDIM, pFin)!=NULL);
            else
                cont = false;

            test_end_file:
            // This is indispensable to correctly handle the end of execution
            // in a "load" command
            while (!cont && act>0 ) {
                fclose(pFin);
                pFin=fstack[--act];
                if (!stopEx && (pFin==stdin)) {
                    interactive=true;
                    cout << ligne << " > ";
                }
                cont = fgets(buffer, BUFDIM, pFin)!=NULL;
            }
            currentPosition=ftell(pFin);
        }
    }
    fclose(pFin);
}

// Insertion d'une commande
void parsefile::insert_command(string command, parsefile *obj, int (*fonc)(parsefile *obj, int argc,char *argv[]))
{
    // Verify we have enough place in the vector

    interpreter p;

    if(command.size()>MAXCMLEN)
        throw parsefile_commandError("Command too long");

    const char *com=command.c_str();
    strncpy(p.commande,com, MAXCMLEN);
    //delete[] com;
    p.fonc=fonc;
    p.obj=obj;

    if(ht_insert(p)==HT_NOTSPACE)
    throw parsefile_commandError("Too many commands");

}

/* Fonction de hashing */
int parsefile::ht_hashing(char com[])
  {
  int cnt=0;
  int somme=0;

  for(cnt=0;com[cnt]!='\0';cnt++)
    somme+=com[cnt];

  somme%=HT_DIM;
  return somme;
}

/* Fonction de recherche du commande dans la table hash  */
int parsefile::ht_find(char com[])
{
  int i,j;
  i=j=ht_hashing(com);
  while(strcmp(com,ht_commande[i].commande)!=0){
    i++;
    i%=HT_DIM;
    if(i==j||ht_commande[i].fonc==NULL) return HT_NOTFOUND;
  }
  return i;
}

/* Fonction d'insertion de la commande */
int parsefile::ht_insert(interpreter newcom)
{
  int i,j;
  i=j=ht_hashing(newcom.commande);
  while(ht_commande[i].fonc!=NULL){
    i++;
    i%=HT_DIM;
    if(i==j) return HT_NOTSPACE;
  }
  ht_commande[i]=newcom;
  return i;
}
