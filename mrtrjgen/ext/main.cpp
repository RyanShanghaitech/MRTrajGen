#include <Python.h>
#include <numpy/arrayobject.h>
#include "Spiral3D.h"

PyObject *GenSpiral3D(PyObject* self, PyObject* args)
{
    bool bRet;
    double dS; double dNp; double dTht0; double dPhi0; double dUTht; double dUPhi; double dDt; double dKmax;
    bRet = PyArg_ParseTuple(args, "dddddddd", &dS, &dNp, &dTht0, &dPhi0, &dUTht, &dUPhi, &dDt, &dKmax);
    if (!bRet) {PyErr_SetString(PyExc_ValueError, "arg parse failed."); return NULL;}

    Spiral3D sp;
    sp.Update(dS, dNp, dTht0, dPhi0, dUTht, dUPhi, dDt, dKmax);
    std::vector<double> vdKx, vdKy, vdKz;
    sp.GetTraj(&vdKx, &vdKy, &vdKz);

    double *data = new double[sp.GetNpt()*3];
    double *ptr = data;
    memcpy(ptr, vdKx.data(), sp.GetNpt()*sizeof(double));
    ptr += sp.GetNpt();
    memcpy(ptr, vdKy.data(), sp.GetNpt()*sizeof(double));
    ptr += sp.GetNpt();
    memcpy(ptr, vdKz.data(), sp.GetNpt()*sizeof(double));
    
    PyArrayObject *pyaTraj;
    {
        npy_intp dims[] = {3, sp.GetNpt()};
        pyaTraj = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, data);
        PyArray_ENABLEFLAGS(pyaTraj, NPY_ARRAY_OWNDATA);
    }
    {
        npy_intp dims[] = {1, 0};
        PyArray_Dims permute = {dims, 2};
        pyaTraj = (PyArrayObject*)PyArray_Transpose(pyaTraj, &permute);
    }

    return (PyObject*)pyaTraj;
}

static PyMethodDef aMeth[] = {
    {"GenSpiral3D", GenSpiral3D, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef sMod = {
    PyModuleDef_HEAD_INIT,
    "mrtrjgen.ext",   /* name of module */
    NULL,
    -1,
    aMeth
};

PyMODINIT_FUNC
PyInit_ext(void) {
    import_array();
    return PyModule_Create(&sMod);
}