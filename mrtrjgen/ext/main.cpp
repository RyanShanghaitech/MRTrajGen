#include <Python.h>
#include <numpy/arrayobject.h>
#include "Spiral3D_A.h"
#include "Spiral3D_B.h"

enum SpiralType
{
    eSpiral3D_A = 0,
    eSpiral3D_B,
};

PyObject *GenSpiral3D(PyObject* self, PyObject* args)
{
    bool bRet;
    int64_t lSpType;
    double dS; double dNp; double dTht0; double dPhi0; double dUTht; double dUPhi; double dDt; double dKmax;
    bRet = PyArg_ParseTuple(args, "Ldddddddd", &lSpType, &dS, &dNp, &dTht0, &dPhi0, &dUTht, &dUPhi, &dDt, &dKmax);
    if (!bRet) {PyErr_SetString(PyExc_ValueError, "arg parse failed."); return NULL;}

    Trajectory *pTraj = NULL;
    if (lSpType == eSpiral3D_A)
    {
        pTraj = new Spiral3D_A;
    }
    else if (lSpType == eSpiral3D_B)
    {
        pTraj = new Spiral3D_B;
    }
    else
    {
        PyErr_SetString(PyExc_ValueError, "Wrong Spiral Type.");
        return NULL;
    }
    pTraj->Update(std::vector<double>({dS, dNp, dTht0, dPhi0, dUTht, dUPhi, dDt, dKmax}));
    std::vector<double> vdKx, vdKy, vdKz;
    pTraj->GetTraj(&vdKx, &vdKy, &vdKz);

    double *data = new double[pTraj->GetNpt()*3];
    double *ptr = data;
    memcpy(ptr, vdKx.data(), pTraj->GetNpt()*sizeof(double));
    ptr += pTraj->GetNpt();
    memcpy(ptr, vdKy.data(), pTraj->GetNpt()*sizeof(double));
    ptr += pTraj->GetNpt();
    memcpy(ptr, vdKz.data(), pTraj->GetNpt()*sizeof(double));
    
    PyArrayObject *pyaTraj;
    {
        npy_intp dims[] = {3, pTraj->GetNpt()};
        pyaTraj = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, data);
        PyArray_ENABLEFLAGS(pyaTraj, NPY_ARRAY_OWNDATA);
    }
    {
        npy_intp dims[] = {1, 0};
        PyArray_Dims permute = {dims, 2};
        pyaTraj = (PyArrayObject*)PyArray_Transpose(pyaTraj, &permute);
    }

    delete pTraj;
    return (PyObject*)pyaTraj;
}

static PyMethodDef aMeth[] = {
    {"GenSpiral3D", GenSpiral3D, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef sMod = {
    PyModuleDef_HEAD_INIT,
    "ext",   /* name of module */
    NULL,
    -1,
    aMeth
};

PyMODINIT_FUNC
PyInit_ext(void) {
    import_array();
    return PyModule_Create(&sMod);
}