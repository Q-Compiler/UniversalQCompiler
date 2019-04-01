#  Copyright 2019 UniversalQCompiler (https://github.com/Q-Compiler/UniversalQCompiler)
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

#include <Python.h>
#include "wstp.h"

// package MathematicaLink

// global variables 
static WSENV mlink_env;
static WSLINK mlink_link;

// helper function to initialize mathematica connection
static int initialize(void) {
    int errn;

    printf("Initializing... ");
    if ((mlink_env = WSInitialize(0))
            == NULL) {
        printf("failed\n");
        // PyErr_SetString()
        return 1;
    }
    printf("OK\n");
    printf("Opening link... ");
    
    if ((mlink_link = WSOpenString(mlink_env, "-linklaunch -linkname math", &errn))
            == NULL) {
        // PyErr_SetString()
        printf("failed\n");
        return 1;
    }
    printf("OK\n");
    return 0;
}

// helper function to import mathematica packages
static int import(const char* str) {
    WSPutFunction(mlink_link, "EvaluatePacket", 1);
     WSPutFunction(mlink_link, "Import", 1);
      WSPutString(mlink_link, str);
    WSEndPacket(mlink_link);

    int pkt;
    while ((pkt = WSNextPacket(mlink_link))) {
        switch (pkt) {
            case MESSAGEPKT:
                printf("Got a warning!\n");
                WSNewPacket(mlink_link);
                return 1;
            case RETURNPKT:
                printf("Got answer\n");
                WSNewPacket(mlink_link);
                return 0;
            default:
                printf("Unknown packet %d\n", pkt);
                WSNewPacket(mlink_link);
        }
    }
    // dummy
    return -1;
}

// function accessed from python
static PyObject *
mlink_import_package(PyObject * self, PyObject *args)
{
    // unpack argument (string!)
    const char * filepath;
    if (!PyArg_ParseTuple(args, "s", &filepath))
        return NULL;

    // import mathematica package
    if (import(filepath) == 1) {
        printf("import failed\n");
        return NULL;
    }
    Py_INCREF(Py_None);
    return Py_None;
}
    
static PyObject *
mlink_evaluate(PyObject * self, PyObject * args)
{
    // unpack argument (string!)
    const char * str;
    if (!PyArg_ParseTuple(args, "s", &str))
        return NULL;

    // Evaluate argument
    WSPutFunction(mlink_link, "EvaluatePacket", 1);
     WSPutFunction(mlink_link, "ToString", 2);
      WSPutFunction(mlink_link, "ToExpression", 1);
       WSPutFunction(mlink_link, "N", 1);
        WSPutString(mlink_link, str);
      WSPutFunction(mlink_link, "ToExpression", 1);
       WSPutString(mlink_link, "InputForm");
    WSEndPacket(mlink_link);
    
    int pkt;
    while ((pkt = (WSNextPacket(mlink_link) != RETURNPKT))) {
        WSNewPacket(mlink_link);
    }
    const char * ans;
    WSGetString(mlink_link, &ans);
    WSNewPacket(mlink_link);

    return PyUnicode_FromString(ans);
}

static PyObject *
mlink_open(PyObject * self, PyObject * args)
{
    // initialize connection to mathematica
    if (initialize() == 1) return NULL;
    
    printf("Activating... ");
    WSActivate(mlink_link);
    printf("OK\n");

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
mlink_close(PyObject * self, PyObject * args)
{
    // closing
    WSClose(mlink_link);
    Py_INCREF(Py_None);
    return Py_None;
}

/*static PyObject *
mlink_atest(PyObject * self, PyObject * args) {*/
    

// necessary Python stuff, irrelevant feature-wise
// details can be found at https://docs.python.org/3.6/extending/extending.html

// method table = functions to export
static PyMethodDef mlink_Methods[] = {
    {"import_package", mlink_import_package, METH_VARARGS, "Import Mathematica packages"},
    {"eval", mlink_evaluate, METH_VARARGS, "Evaluate mathematica expression"},
    {"open", mlink_open,   METH_VARARGS, "Open Mathematica link"},
    {"close", mlink_close,   METH_VARARGS, "Close Mathematica link"},
    {NULL, NULL, 0, NULL}
};

// module definition
static struct PyModuleDef mlinkmodule = {
    PyModuleDef_HEAD_INIT,
    "MathematicaLink",
    NULL,
    -1,
    mlink_Methods
};

// initialisation function (only non-static)
PyMODINIT_FUNC PyInit_MathematicaLink(void) {

    return PyModule_Create(&mlinkmodule);
}
