/***************************************************************************
*     CLASS parsefile                                                      *
*     Davide Bucci, CROMA     2005-2012                                    *
*     MINATEC-INPG, 3, parvis Luis Neel                                    *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.inpg.fr                                                *
*                                                                          *
*     Version: 1.3                                                         *
****************************************************************************
*    C++ file parsing engine                                               *
****************************************************************************/

/* Version 1.3, January 2012 */
/* Version 1.2, January 2009 */
/* Version 1.1, February 2006 */
/* Version 1.0, 2005 */


#ifndef PARSEFILE_H
#define PARSEFILE_H

#include <string>
#include <cstring>
#include <iostream>
#include <map>
#include "parsercl.h"
#include <stack>

using namespace std;

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

#define HT_NOTFOUND -1
#define HT_NOTSPACE -1
#define MAXCMLEN     64             // Longueur maximale d'un nom de commande
#define BUFDIM      256             // Longueur maximale d'une ligne du fichier
#define MAXFILE     20              // Nombre max de files ouverts

typedef class parsefile parsefile_tag;

typedef struct tag_interpreter{
      char commande[MAXCMLEN];  /* clé */
      parsefile_tag *obj;           // Objet sur lequel executer les fonctions
                                    // (qui doivent être statiques)
      int (* fonc)(parsefile_tag *obj, int argc,char *argv[]); /* fonction à
                                    laquelle sauter */
} interpreter;


class parsefile {

    /* Dimension de l'hash table pour stocker les commandes à executer */
    #define HT_DIM 100

    /* Hash table containing labels */
    map<string, long int> labels;
    map<string, double> for_end;
    map<string, double> for_next;
    map<string, double> for_increment;

    bool noForCycles;

    /* Hash table contenant l'adresse des fonctions auxquelles sauter pour executer la commande */
    interpreter ht_commande[HT_DIM];

    virtual int ht_hashing(char com[]);
    virtual int ht_find(char com[]);
    virtual int ht_insert(interpreter newcom);

    FILE *fstack[MAXFILE];
    int act;
    FILE *pFin;
    long int currentPosition;
    long int successivePosition;

    bool system_command_allowed;

    // Keep track of the last read in a file. True if the last read
    // has been successful, false otherwise.
    bool lastRead;

    // Add some spaces between arguments of print and fprint
    // as well as newlines.
    bool addspace;

    // Progressive number of 'if' instructions in the file.
    stack<string> ifstack;
    // Progressive number of 'for' instructions in the file.
    stack<string> forstack;

    // List of opened files.
    map<double, FILE*> openedfiles;
    double filecounter;

protected:
    // Functions to get and set a bookmark.
    long int getBookmark(void);
    void setBookmark(long int b);
    numParser *nP;
public:

    static int c_load(parsefile *obj, int argc,char *argv[]);
    static int c_system(parsefile *obj, int argc,char *argv[]);

    static int c_goto(parsefile *obj, int argc,char *argv[]);
    static int c_label(parsefile *obj, int argc,char *argv[]);

    static int c_for(parsefile *obj, int argc,char *argv[]);
    static int c_next(parsefile *obj, int argc,char *argv[]);
    static int c_do(parsefile *obj, int argc,char *argv[]);
    static int c_while(parsefile *obj, int argc,char *argv[]);

    static int c_if(parsefile *obj, int argc,char *argv[]);
    static int c_else(parsefile *obj, int argc,char *argv[]);
    static int c_endif(parsefile *obj, int argc,char *argv[]);

    static int c_addspace(parsefile *obj, int argc,char *argv[]);

    static int c_print(parsefile *obj, int argc,char *argv[]);

    static int c_fopen(parsefile *obj, int argc,char *argv[]);
    static int c_fclose(parsefile *obj, int argc,char *argv[]);
    static int c_fprint(parsefile *obj, int argc,char *argv[]);
    static int c_fscan(parsefile *obj, int argc,char *argv[]);
    static int c_fcheck(parsefile *obj, int argc,char *argv[]);
    static int c_fskip(parsefile *obj, int argc,char *argv[]);

    // Lecture d'un fichier
    virtual void read(string fileName, bool interactive);
    // Lecture d'un fichier (à partir d'un pointeur à FILE)
    void readFile(FILE* pFin, bool interactive);


    // Insertion d'une commande
    virtual void insert_command(string command, parsefile_tag *obj,
        int (*fonc)(parsefile *obj, int argc,char *argv[]));

    // Constructeur standard
    parsefile(numParser *p);

    void allow_system_command(bool r)
    {
        system_command_allowed=true;
    }

    void insertVar(string name, double d)
    {
        nP->insertVar(name.c_str(), d);
    }
};


// Structure lancée comme exception si une commande n'a pas eu succés
class parsefile_commandError:exception {
//  static const int MAXLEN = 100;
//  char errMess[MAXLEN];
    string errMess;
public:

    parsefile_commandError()
    {
        errMess="Error";
        //strncat(errMess, "Error\n",MAXLEN);
    }
    parsefile_commandError(string m)
    {
        //cout<<m; cout.flush();
        /*errMess[m.getLen()]='\0';
        // HEY THIS CODE DOES NOT ADD \0!!!
        strncat(errMess, m.c_str(),MAXLEN);*/
        cerr << m<<endl;
        errMess=m;
    }
    ~parsefile_commandError() throw () {};

    const char* getMess()
    {
        return errMess.c_str();
    }

};

class parsefile_stop:exception {};
#endif
