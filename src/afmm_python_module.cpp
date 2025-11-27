
/***************************************************************************
*     RCWA/FMM 3D full vectorial mode solver and propagator                *
*     Pyton interface                                                      *
*                                                                          *
*     Davide Bucci, CROMA                                                  *
*     Jérôme Michallon, CROMA                                              *
*     MINATEC-INPG, 3, parvis Louis Neel                                   *
*     38016, Grenoble CEDEX, France                                        *
*                                                                          *
*     bucci@minatec.grenoble-inp.fr                                        *
*                                                                          *
****************************************************************************/

// g++ -dynamiclib -I/usr/include/python/ -lpython -o pyAFMM.so afmm_python_module.cpp

#if PY_MAJOR_VERSION >= 3
    #error Python 3 is requested for the compilation
#endif

#include <Python.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>

using namespace std;

#include "structure.h"
#include "block_matrix.h"
#include "commands.h"
#include "compileinfo.h"

extern "C" {
    PyMODINIT_FUNC PyInit_pyAFMM(void);
    static PyObject* py_AFMM_banner(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_parsescript(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_size(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_harmonics(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_wavelength(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_substrate(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_solve(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_rectangle(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_clear(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_lowindex(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_highindex(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_carpet(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_assemble(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_pml(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_pml_transf(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_bend(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_inpstruct(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_outgmodes(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_indmatrix(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_wants(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_section(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_select(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_order(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_bloch(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_spectrum(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_powerz(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_monitor(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_angles(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_coefficient(PyObject* self, PyObject* args);
    static PyObject* py_AFMM_coefficient_id(PyObject* self, PyObject* args);
}

/*
 * Bind Python function names to our C functions
 */
static PyMethodDef myModule_methods[] = {
    {"banner", py_AFMM_banner, METH_VARARGS},
    {"parsescript", py_AFMM_parsescript, METH_VARARGS},
    {"size", py_AFMM_size, METH_VARARGS},
    {"harmonics", py_AFMM_harmonics, METH_VARARGS},
    {"wavelength", py_AFMM_wavelength, METH_VARARGS},
    {"substrate", py_AFMM_substrate, METH_VARARGS},
    {"solve", py_AFMM_solve, METH_VARARGS},
    {"rectangle", py_AFMM_rectangle, METH_VARARGS},
    {"clear", py_AFMM_clear, METH_VARARGS},
    {"lowindex", py_AFMM_lowindex, METH_VARARGS},
    {"highindex", py_AFMM_highindex, METH_VARARGS},
    {"carpet", py_AFMM_carpet, METH_VARARGS},
    {"assemble", py_AFMM_assemble, METH_VARARGS},
    {"pml", py_AFMM_pml, METH_VARARGS},
    {"pml_transf", py_AFMM_pml_transf, METH_VARARGS},
    {"bend", py_AFMM_bend, METH_VARARGS},
    {"inpstruct", py_AFMM_inpstruct, METH_VARARGS},
    {"outgmodes", py_AFMM_outgmodes, METH_VARARGS},
    {"indmatrix", py_AFMM_indmatrix, METH_VARARGS},
    {"wants", py_AFMM_wants, METH_VARARGS},
    {"section", py_AFMM_section, METH_VARARGS},
    {"select", py_AFMM_select, METH_VARARGS},
    {"order", py_AFMM_order, METH_VARARGS},
    {"bloch", py_AFMM_bloch, METH_VARARGS},
    {"spectrum", py_AFMM_spectrum, METH_VARARGS},
    {"powerz", py_AFMM_powerz, METH_VARARGS},
    {"monitor", py_AFMM_monitor, METH_VARARGS},
    {"angles", py_AFMM_angles, METH_VARARGS},
    {"coefficient", py_AFMM_coefficient, METH_VARARGS},
    {"coefficient_id", py_AFMM_coefficient_id, METH_VARARGS},

    {NULL, NULL,0,NULL}
};

static PyObject *afmmError;

static struct PyModuleDef pyAFMMmodule = {
    PyModuleDef_HEAD_INIT,
    "pyAFMM",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    myModule_methods
};

numParser nP;

structure waveguide(&nP);
commands co(waveguide, &nP);
bool externalCommands = false;

/*
 * Python2 calls this to let us initialize our module
 *
void initpyAFMM(void)
{
    PyObject *m = Py_InitModule("pyAFMM", myModule_methods);
    if (m == NULL)
        return;

    afmmError = PyErr_NewException("pyAFMM.error", NULL, NULL);
    Py_INCREF(afmmError);
    PyModule_AddObject(m, "error", afmmError);
}
*/

/*
 * Python3 calls this to let us initialize our module
 */
PyMODINIT_FUNC PyInit_pyAFMM(void)
{
    PyObject *m;
    CHECK_EXPIRED();

    m = PyModule_Create(&pyAFMMmodule);
    if (m == NULL)
        return NULL;

    afmmError = PyErr_NewException("pyAFMM.error", NULL, NULL);
    Py_XINCREF(afmmError);
    if (PyModule_AddObject(m, "error", afmmError) < 0) {
        Py_XDECREF(afmmError);
        Py_CLEAR(afmmError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}


/**
  Parse an AFMM script.
  Execution of external commands is active.
*/
static PyObject* py_AFMM_parsescript(PyObject* self, PyObject* args)
{
    const char *script;

    if (!PyArg_ParseTuple(args, "s", &script))
        return NULL;

    co.allow_system_command(true);
    FILE *f= tmpfile();
    fprintf(f, "%s",script);

    try {
        co.readFile(f,false);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    //fclose(f);
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  size command
*/
static PyObject* py_AFMM_size(PyObject* self, PyObject* args)
{
    double sx, sy;

    if (!PyArg_ParseTuple(args, "dd", &sx, &sy))
        return NULL;

    try {
        waveguide.do_size(sx,sy);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  harmonics command
*/
static PyObject* py_AFMM_harmonics(PyObject* self, PyObject* args)
{
    int hx, hy;

    if (!PyArg_ParseTuple(args, "ii", &hx, &hy))
        return NULL;

    try {
        waveguide.do_harmonics(hx,hy);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}


/**
  wavelength command
*/
static PyObject* py_AFMM_wavelength(PyObject* self, PyObject* args)
{
    double l;

    if (!PyArg_ParseTuple(args, "d", &l))
        return NULL;

    waveguide.set_wavelength(l);

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  substrate command
*/
static PyObject* py_AFMM_substrate(PyObject* self, PyObject* args)
{
    double r,i;
    Py_complex c;

    if (!PyArg_ParseTuple(args, "D", &c))
        return NULL;

    r=c.real;
    i=c.imag;

    waveguide.getCur()->set_substrate(complex<double>(r,i));

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  solve command
  Return: a Python list of complex containing the effective indices of the
    last section in the structure.
*/
static PyObject* py_AFMM_solve(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    vector<double> results;
    try {
        results = waveguide.do_solve();
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    // results vector contains complex data stacked in real part followed by
    // imaginary part. Hence, the number of elements should be divided by 2.

    PyObject *returnObj = PyList_New(results.size()/2);
    if (returnObj==NULL) {
        PyErr_SetString(afmmError, "Could not create a Python list from C++.");
        return NULL;
    }

    unsigned int k=0;
    for(unsigned int i=0;i<results.size();i+=2) {
        PyObject *cdata = PyComplex_FromDoubles(results[i], results[i+1]);
        if (cdata==NULL) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not create a complex from C++.");
            return NULL;
        }
        if(PyList_SetItem(returnObj, k++, cdata)) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not add to list from C++.");
            return NULL;
        }
    }
    return returnObj;
}

/**
  Show the program banner
 */
static PyObject* py_AFMM_banner(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    cout << " ***************************************************************************\n"
         << " *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *\n"
         << " *                            version 1.4.6                                *\n"
         << " *                                                                         *\n"
         << " *     Build date: " << __DATE__<<  "                                             *\n"
         << " *     Source revision: "
         << setw(30)
         << setfill(' ')
         << left
         << _SVNGLOBALVERSION
         << std::setw(0)
         << "                     *\n"
         << " *                                                                         *\n"
         << " *     Davide Bucci, CROMA     March 2008 - current                        *\n"
         << " *     Jérôme Michallon, CROMA     May 2012 - February 2014                *\n"
         << " *     MINATEC-Grenoble INP, 3, parvis Louis Neel                          *\n"
         << " *     38016, Grenoble CEDEX, France                                       *\n"
         << " *                                                                         *\n"
         << " *     bucci@minatec.grenoble-inp.fr                                       *\n"
         << " *                                                                         *\n"
         << " ***************************************************************************\n";

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  rectangle command
*/
static PyObject* py_AFMM_rectangle(PyObject* self, PyObject* args)
{
    double r,i;
    Py_complex c;
    double wx, wy, px, py;

    if (!PyArg_ParseTuple(args, "Ddddd", &c,&wx,&wy,&px,&py))
        return NULL;

    r=c.real;
    i=c.imag;

    waveguide.getCur()->add_rectangle(complex<double>(r,i), wx, wy, px, py);

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  lowindex command
*/
static PyObject* py_AFMM_lowindex(PyObject* self, PyObject* args)
{
    Py_complex c;

    if (!PyArg_ParseTuple(args, "D", &c))
        return NULL;

    waveguide.getCur()->set_lowindex(complex<double>(c.real,c.imag));

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  highindex command
*/
static PyObject* py_AFMM_highindex(PyObject* self, PyObject* args)
{
    Py_complex c;

    if (!PyArg_ParseTuple(args, "D", &c))
        return NULL;
    waveguide.getCur()->set_highindex(complex<double>(c.real,c.imag));

    Py_INCREF(Py_None);
    return Py_None;
}


/**
  clear command
*/
static PyObject* py_AFMM_clear(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    waveguide.reset();
    cout<<"The structure definition has been cleared.\n";
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  carpet command
*/

static PyObject* py_AFMM_carpet(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    waveguide.set_ensureConvergence(true);

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  assemble command
*/

static PyObject* py_AFMM_assemble(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;
    try {
        waveguide.do_assemble();
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  pml command
*/
static PyObject* py_AFMM_pml(PyObject* self, PyObject* args)
{
    double wx,wy;
    Py_complex alpha;

    if (!PyArg_ParseTuple(args, "Ddd", &alpha, &wx,&wy))
        return NULL;

    try {
        waveguide.getCur()->set_pml(wx,wy,
            complex<double>(alpha.real, alpha.imag));
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  pml_trasnf command
*/
static PyObject* py_AFMM_pml_transf(PyObject* self, PyObject* args)
{
    double qdx,qdy;
    Py_complex g;

    if (!PyArg_ParseTuple(args, "ddD", &qdx,&qdy,&g))
        return NULL;

    double r=g.real;
    double i=g.imag;

    try {
        waveguide.getCur()->set_pml_transf(qdx,qdy,complex<double>(r,i));
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  bend command
*/
static PyObject* py_AFMM_bend(PyObject* self, PyObject* args)
{
    double r;

    if (!PyArg_ParseTuple(args, "d", &r))
        return NULL;

    try {
        waveguide.getCur()->set_bend(r);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  inpstruct command
*/
static PyObject* py_AFMM_inpstruct(PyObject* self, PyObject* args)
{
    double sx, sy;
    char *otype;
    db_matrix out;

    if (!PyArg_ParseTuple(args, "ddz", &sx, &sy, &otype))
        return NULL;

    try {
        out = waveguide.getCur()->do_inpstruct(sx, sy, otype);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    PyObject *returnObj = PyList_New(sy);
    for(unsigned int i=0; i<sy;++i) {
        PyObject *rowObj = PyList_New(sx);
        if (returnObj==NULL) {
            PyErr_SetString(afmmError,
                "Could not create a Python list from C++.");
            return NULL;
        }

        for(unsigned int j=0;j<sx;++j) {
            PyObject *cdata =
                PyComplex_FromDoubles(out(i,j).real(),out(i,j).imag());
            if (cdata==NULL) {
                Py_DECREF(returnObj);
                PyErr_SetString(afmmError,
                    "Could not create a complex from C++.");
                return NULL;
            }
            if(PyList_SetItem(rowObj, j, cdata)) {
                Py_DECREF(rowObj);
                PyErr_SetString(afmmError, "Could not add to list from C++.");
                return NULL;
            }
        }
        if(PyList_SetItem(returnObj, i, rowObj)) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not add to list from C++.");
            return NULL;
        }
    }
    return returnObj;
}

/**
  outgmodes command
*/
static PyObject* py_AFMM_outgmodes(PyObject* self, PyObject* args)
{
    double sx, sy;
    char *fieldcomponent;
    list<db_matrix> modes;

    if (!PyArg_ParseTuple(args, "zdd", &fieldcomponent, &sx, &sy))
        return NULL;
    
    modes =  waveguide.getCur()->do_outgmodes(sx, sy, fieldcomponent, S, "");
    printf("Matrix calculated correctly\n");
    fflush(stdout);

    PyObject *returnObj = PyList_New(modes.size());
    int k=0;
    for (list<db_matrix>::iterator out = modes.begin(); out!=modes.end();++out)
    {
        PyObject *colObj = PyList_New(sy);
        for(unsigned int i=0; i<sy;++i) {
            PyObject *rowObj = PyList_New(sx);
            for(unsigned int j=0;j<sx;++j) {
                PyObject *cdata =
                    PyComplex_FromDoubles((*out)(i,j).real(),
                        (*out)(i,j).imag());

                if(PyList_SetItem(rowObj, j, cdata)) {
                    Py_DECREF(rowObj);
                    PyErr_SetString(afmmError,
                        "Could not add to list from C++.");
                    return NULL;
                }
            }
            if(PyList_SetItem(colObj, i, rowObj)) {
                Py_DECREF(colObj);
                PyErr_SetString(afmmError, "Could not add to list from C++.");
                return NULL;
            }
        }
        if(PyList_SetItem(returnObj, k++, colObj)) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not add to list from C++.");
            return NULL;
        }
    }
    return returnObj;
}

/**
  indmatrix command. The Python-esque equivalent of the indfile matrix that
  returns a matrix and do not write on a file.
*/
static PyObject* py_AFMM_indmatrix(PyObject* self, PyObject* args)
{
    PyObject *rows;
    if (!PyArg_ParseTuple(args, "O", &rows)) {
        PyErr_SetString(afmmError, "Can not obtain an object");
        return NULL;
    }
    if(!PyList_Check(rows)) {
        PyErr_SetString(afmmError, "Object is not a list");
        return NULL;
    }
    if(!PyList_Check(PyList_GetItem(rows,0))) {
        PyErr_SetString(afmmError, "Object is not a list of lists");
        return NULL;
    }

    unsigned int nrow = PyList_Size(rows);
    unsigned int ncol = PyList_Size(PyList_GetItem(rows,0));

    db_matrix rd(nrow, ncol);

    PyObject *col;
    PyObject *v;

    for(unsigned int i=0; i<nrow; ++i){
        col=PyList_GetItem(rows,i);
        for(unsigned int j=0; j<ncol; ++j) {
            v=PyList_GetItem(col,j);
            if(!PyComplex_Check(v)) {
                PyErr_SetString(afmmError, "Elements are not complex");
                return NULL;
            }
            rd(i,j)=complex<double>(PyComplex_RealAsDouble(v),
                PyComplex_ImagAsDouble(v));
        }
    }
    try {
        waveguide.getCur()->store_refractive_index(rd);
        cout<<"Refractive index matrix loaded."<<endl;
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/**
  wants command
*/
static PyObject* py_AFMM_wants(PyObject* self, PyObject* args)
{
    const char *code;

    if (!PyArg_ParseTuple(args, "s", &code))
        return NULL;
    try {
        waveguide.do_wants(code);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}


/**
  section command
*/
static PyObject* py_AFMM_section(PyObject* self, PyObject* args)
{
    double tl;

    if (!PyArg_ParseTuple(args, "d", &tl))
        return NULL;

    try {
        waveguide.do_section(tl);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}


/**
  select command
*/
static PyObject* py_AFMM_select(PyObject* self, PyObject* args)
{
    unsigned int number;

    if (!PyArg_ParseTuple(args, "I", &number))
        return NULL;
    // The internal numbering starts from 0, as with C/C++

     try {
        waveguide.do_select(number);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  order command
*/
static PyObject* py_AFMM_order(PyObject* self, PyObject* args)
{
    double minv, maxv;

    if (!PyArg_ParseTuple(args, "dd", &minv, &maxv))
        return NULL;

    try {
        waveguide.getCur()->do_order(minv,maxv);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  bloch command
  Return: a Python list of complex containing the effective indices
*/
static PyObject* py_AFMM_bloch(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    vector<double> results;
    try {
        results = waveguide.do_bloch();
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    // results vector contains complex data stacked in real part followed by
    // imaginary part. Hence, the number of elements should be divided by 2.

    PyObject *returnObj = PyList_New(results.size()/2);
    if (returnObj==NULL) {
        PyErr_SetString(afmmError, "Could not create a Python list from C++.");
        return NULL;
    }

    unsigned int k=0;
    for(unsigned int i=0;i<results.size();i+=2) {
        PyObject *cdata = PyComplex_FromDoubles(results[i], results[i+1]);
        if (cdata==NULL) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not create a complex from C++.");
            return NULL;
        }
        if(PyList_SetItem(returnObj, k++, cdata)) {
            Py_DECREF(returnObj);
            PyErr_SetString(afmmError, "Could not add to list from C++.");
            return NULL;
        }
    }
    return returnObj;
}

/**
  spectrum command
  Return: a Python list of complex containing the effective indices of the
    last section in the structure.
*/
static PyObject* py_AFMM_spectrum(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    PyObject *returnObj;

    try {
        // Check if the C matrix is empty
        if (waveguide.getCur()->getW().isEmpty()) {
            throw parsefile_commandError(
                "spectrum: the problem must be set and modes calculated with"
                "solve before trying to write the spectrum on a file.");
        }
        returnObj = PyList_New(waveguide.getCur()->getB().getNrow());
        if (returnObj==NULL) {
            PyErr_SetString(afmmError,
                "Could not create a Python list from C++.");
            return NULL;
        }
        int k=0;
        for(int i=0; i<waveguide.getCur()->getB().getNrow(); ++i) {
            complex<double>e=waveguide.getCur()->getB()(i,i);
            e = waveguide.getEffectiveIndex(e);
            PyObject *cdata = PyComplex_FromDoubles(e.real(), e.imag());
            if (cdata==NULL) {
                Py_DECREF(returnObj);
                PyErr_SetString(afmmError,
                    "Could not create a complex from C++.");
                return NULL;
            }
            if(PyList_SetItem(returnObj, k++, cdata)) {
                Py_DECREF(returnObj);
                PyErr_SetString(afmmError,
                    "Could not add to list from C++.");
                return NULL;
            }
        }
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    return returnObj;
}

/**
  powerz command
  Return: a Python floating point containing the calculated power.
*/
static PyObject* py_AFMM_powerz(PyObject* self, PyObject* args)
{
    double z=0;
    double power;

    if (!PyArg_ParseTuple(args, "d", &z))
        return NULL;

    try {
        power = waveguide.do_powerz(z);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    // results vector contains complex data stacked in real part followed by
    // imaginary part. Hence, the number of elements should be divided by 2.

    
    return Py_BuildValue("d", power);
}

/**
  monitor command
  Return: a Python floating point containing the calculated power.
*/
static PyObject* py_AFMM_monitor(PyObject* self, PyObject* args)
{
    double z, wx, wy, px, py;
    double power;

    if (!PyArg_ParseTuple(args, "ddddd", &z, &wx, &wy, &px, &py))
        return NULL;

    try {
        power = waveguide.do_monitor(z,wx,wy,px,py);
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    // results vector contains complex data stacked in real part followed by
    // imaginary part. Hence, the number of elements should be divided by 2.

    return Py_BuildValue("d", power);
}

/**
  angles command
*/
static PyObject* py_AFMM_angles(PyObject* self, PyObject* args)
{
    double n0,thetax, thetay;

    if (!PyArg_ParseTuple(args, "ddd", &n0, &thetax, &thetay))
        return NULL;

    try {
        waveguide.setAngles(n0, thetax, thetay);
        cout << "Angles set to: "<<thetax<<" rad and "<<thetay<<
            " rad in a section with refractive index "<<n0<<"\n";
    } catch (parsefile_commandError e) {
        PyErr_SetString(afmmError, e.getMess());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

/**
  coefficient command
*/
static PyObject* py_AFMM_coefficient(PyObject* self, PyObject* args)
{
    Py_complex coeff;
    char *type;
    char *second;
    double value;
    int iMode;
    complex<double> e;

    if(!PyArg_ParseTuple(args, "zd", &type, &value)) {
        PyErr_SetString(afmmError, "coefficient: could not read the "
                "parameters.\n");
            return NULL;
    }

   iMode = waveguide.getCur()->select_real(value);

    if (iMode<0) {
        cout<<"I could not find the given mode.\n";
        return NULL;
    }
    cout << "Selected mode: ";
    cout << waveguide.getEffectiveIndex(waveguide.getCur()->B(iMode,iMode)) 
        << "\n";

    if(strcmp(type,"f")==0) {
        if(waveguide.getCur()->sWp.isEmpty()) {
            PyErr_SetString(afmmError, "coefficient: excitation coefficients"
                " have not been calculated yet.\n");
            return NULL;
        } else {
            cout<<"Forward ";
            e = waveguide.getCur()->sWp(iMode, 0);
        }
    } else if(strcmp(type,"b")==0) {
        if(waveguide.getCur()->sWm.isEmpty()) {
            PyErr_SetString(afmmError, "coefficient: excitation coefficients"
                " have not been calculated yet.\n");
            return NULL;
        } else {
            cout<<"Backward ";
            e = waveguide.getCur()->sWm(iMode, 0);
        }
    } else {
        PyErr_SetString(afmmError, "coefficient: unrecognized direction "
            "of the excitation. It must be {b|f}.\n");
        return NULL;
    }
    cout << "mode coefficient: " << e << "\n";
    cout << "Squared module: " << abs(e)*abs(e) << "\n";

    // Save the calculated data as an array in the output variable "ans".
    double *s =new double[3];
    coeff.real=e.real();
    coeff.imag=e.imag();


    return Py_BuildValue("D", &coeff);
}

/**
  coefficient_id command
*/
static PyObject* py_AFMM_coefficient_id(PyObject* self, PyObject* args)
{
    Py_complex coeff;
    char *type;
    char *second;
    double value;
    int iMode;
    complex<double> e;

    if(!PyArg_ParseTuple(args, "zi", &type, &iMode)) {
        PyErr_SetString(afmmError, "coefficient_id: could not read the "
                "parameters.\n");
            return NULL;
    }

    if (iMode<0) {
        cout<<"Mode index invalid.\n";
        return NULL;
    }
    cout << "Selected mode: ";
    cout << waveguide.getEffectiveIndex(waveguide.getCur()->B(iMode,iMode)) 
        << "\n";

    if(strcmp(type,"f")==0) {
        if(waveguide.getCur()->sWp.isEmpty()) {
            PyErr_SetString(afmmError, "coefficient_id: excitation coefficients"
                " have not been calculated yet.\n");
            return NULL;
        } else {
            cout<<"Forward ";
            e = waveguide.getCur()->sWp(iMode, 0);
        }
    } else if(strcmp(type,"b")==0) {
        if(waveguide.getCur()->sWm.isEmpty()) {
            PyErr_SetString(afmmError, "coefficient_id: excitation coefficients"
                " have not been calculated yet.\n");
            return NULL;
        } else {
            cout<<"Backward ";
            e = waveguide.getCur()->sWm(iMode, 0);
        }
    } else {
        PyErr_SetString(afmmError, "coefficient_id: unrecognized direction "
            "of the excitation. It must be {b|f}.\n");
        return NULL;
    }
    cout << "mode coefficient: " << e << "\n";
    cout << "Squared module: " << abs(e)*abs(e) << "\n";

    // Save the calculated data as an array in the output variable "ans".
    double *s =new double[3];
    coeff.real=e.real();
    coeff.imag=e.imag();


    return Py_BuildValue("D", &coeff);
}

