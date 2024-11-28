#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include "Spiral2D.h"
#include "Spiral3D_A.h"
#include "Spiral3D_B.h"

PyObject *GenSpiral2D(PyObject* self, PyObject* args)
{
    bool bRet;
    double dKtht1, dKtht2, dTht0, dKmax, dS, dDt, dOv;
    bRet = PyArg_ParseTuple(args, "ddddddd", &dKtht1, &dKtht2, &dTht0, &dKmax, &dS, &dDt, &dOv);
    if (!bRet) {PyErr_SetString(PyExc_ValueError, "arg parse failed."); return NULL;}

    // construct and update trajectory object
    Spiral2D spiral;
    spiral.Update(dKtht1, dKtht2, dTht0, dKmax, dS, dDt, dOv);

    // copy value from spiral object
    std::vector<double> vdKx, vdKy;
    spiral.GetTraj(&vdKx, &vdKy, NULL);
    std::vector<float> vfGx, vfGy;
    spiral.GetGrad(&vfGx, &vfGy, NULL);

    double *adTraj = new double[spiral.GetNpt()*2];
    float *adGrad = new float[spiral.GetNpt()*2];
    {
        double *ptr = adTraj;
        memcpy(ptr, vdKx.data(), spiral.GetNpt()*sizeof(double));
        ptr += spiral.GetNpt();
        memcpy(ptr, vdKy.data(), spiral.GetNpt()*sizeof(double));
    }
    {
        float *ptr = adGrad;
        memcpy(ptr, vfGx.data(), spiral.GetNpt()*sizeof(float));
        ptr += spiral.GetNpt();
        memcpy(ptr, vfGy.data(), spiral.GetNpt()*sizeof(float));
    }
    
    // construct PyArrayObject
    PyArrayObject *pyaTraj, *pyaGrad;
    {
        npy_intp dims[] = {2, spiral.GetNpt()};
        pyaTraj = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, adTraj);
        PyArray_ENABLEFLAGS(pyaTraj, NPY_ARRAY_OWNDATA);
        pyaGrad = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, adGrad);
        PyArray_ENABLEFLAGS(pyaGrad, NPY_ARRAY_OWNDATA);
    }
    {
        npy_intp dims[] = {1, 0};
        PyArray_Dims permute = {dims, 2};
        pyaTraj = (PyArrayObject*)PyArray_Transpose(pyaTraj, &permute);
        pyaGrad = (PyArrayObject*)PyArray_Transpose(pyaGrad, &permute);
    }

    return Py_BuildValue("OO", pyaTraj, pyaGrad);
}

PyObject *GenSpiral3D_A(PyObject* self, PyObject* args)
{
    bool bRet;
    double dNp, dUTht, dUPhi, dTht0, dPhi0, dKmax, dS, dDt;
    bRet = PyArg_ParseTuple(args, "dddddddd", &dNp, &dUTht, &dUPhi, &dTht0, &dPhi0, &dKmax, &dS, &dDt);
    if (!bRet) {PyErr_SetString(PyExc_ValueError, "arg parse failed."); return NULL;}

    // construct and update trajectory object
    Spiral3D_A spiral;
    spiral.Update(dNp, dUTht, dUPhi, dTht0, dPhi0, dKmax, dS, dDt);

    // copy value from spiral object
    std::vector<double> vdKx, vdKy, vdKz;
    spiral.GetTraj(&vdKx, &vdKy, &vdKz);
    std::vector<float> vfGx, vfGy, vfGz;
    spiral.GetGrad(&vfGx, &vfGy, &vfGz);

    double *adTraj = new double[spiral.GetNpt()*3];
    float *adGrad = new float[spiral.GetNpt()*3];
    {
        double *ptr = adTraj;
        memcpy(ptr, vdKx.data(), spiral.GetNpt()*sizeof(double));
        ptr += spiral.GetNpt();
        memcpy(ptr, vdKy.data(), spiral.GetNpt()*sizeof(double));
        ptr += spiral.GetNpt();
        memcpy(ptr, vdKz.data(), spiral.GetNpt()*sizeof(double));
    }
    {
        float *ptr = adGrad;
        memcpy(ptr, vfGx.data(), spiral.GetNpt()*sizeof(float));
        ptr += spiral.GetNpt();
        memcpy(ptr, vfGy.data(), spiral.GetNpt()*sizeof(float));
        ptr += spiral.GetNpt();
        memcpy(ptr, vfGz.data(), spiral.GetNpt()*sizeof(float));
    }
    
    // construct PyArrayObject
    PyArrayObject *pyaTraj, *pyaGrad;
    {
        npy_intp dims[] = {3, spiral.GetNpt()};
        pyaTraj = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, adTraj);
        PyArray_ENABLEFLAGS(pyaTraj, NPY_ARRAY_OWNDATA);
        pyaGrad = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, adGrad);
        PyArray_ENABLEFLAGS(pyaGrad, NPY_ARRAY_OWNDATA);
    }
    {
        npy_intp dims[] = {1, 0};
        PyArray_Dims permute = {dims, 2};
        pyaTraj = (PyArrayObject*)PyArray_Transpose(pyaTraj, &permute);
        pyaGrad = (PyArrayObject*)PyArray_Transpose(pyaGrad, &permute);
    }

    return Py_BuildValue("OO", pyaTraj, pyaGrad);
}

PyObject *GenSpiral3D_B(PyObject* self, PyObject* args)
{
    bool bRet;
    double dNp, dUTht, dUPhi, dTht0, dPhi0, dKmax, dS, dDt;
    bRet = PyArg_ParseTuple(args, "dddddddd", &dNp, &dUTht, &dUPhi, &dTht0, &dPhi0, &dKmax, &dS, &dDt);
    if (!bRet) {PyErr_SetString(PyExc_ValueError, "arg parse failed."); return NULL;}

    // construct and update trajectory object
    Spiral3D_B spiral;
    spiral.Update(dNp, dUTht, dUPhi, dTht0, dPhi0, dKmax, dS, dDt);

    // copy value from spiral object
    std::vector<double> vdKx, vdKy, vdKz;
    spiral.GetTraj(&vdKx, &vdKy, &vdKz);
    std::vector<float> vfGx, vfGy, vfGz;
    spiral.GetGrad(&vfGx, &vfGy, &vfGz);

    double *adTraj = new double[spiral.GetNpt()*3];
    float *adGrad = new float[spiral.GetNpt()*3];
    {
        double *ptr = adTraj;
        memcpy(ptr, vdKx.data(), spiral.GetNpt()*sizeof(double));
        ptr += spiral.GetNpt();
        memcpy(ptr, vdKy.data(), spiral.GetNpt()*sizeof(double));
        ptr += spiral.GetNpt();
        memcpy(ptr, vdKz.data(), spiral.GetNpt()*sizeof(double));
    }
    {
        float *ptr = adGrad;
        memcpy(ptr, vfGx.data(), spiral.GetNpt()*sizeof(float));
        ptr += spiral.GetNpt();
        memcpy(ptr, vfGy.data(), spiral.GetNpt()*sizeof(float));
        ptr += spiral.GetNpt();
        memcpy(ptr, vfGz.data(), spiral.GetNpt()*sizeof(float));
    }
    
    // construct PyArrayObject
    PyArrayObject *pyaTraj, *pyaGrad;
    {
        npy_intp dims[] = {3, spiral.GetNpt()};
        pyaTraj = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT64, adTraj);
        PyArray_ENABLEFLAGS(pyaTraj, NPY_ARRAY_OWNDATA);
        pyaGrad = (PyArrayObject*)PyArray_SimpleNewFromData(2, dims, NPY_FLOAT32, adGrad);
        PyArray_ENABLEFLAGS(pyaGrad, NPY_ARRAY_OWNDATA);
    }
    {
        npy_intp dims[] = {1, 0};
        PyArray_Dims permute = {dims, 2};
        pyaTraj = (PyArrayObject*)PyArray_Transpose(pyaTraj, &permute);
        pyaGrad = (PyArrayObject*)PyArray_Transpose(pyaGrad, &permute);
    }

    return Py_BuildValue("OO", pyaTraj, pyaGrad);
}

static PyMethodDef aMeth[] = {
    {"GenSpiral2D", GenSpiral2D, METH_VARARGS, ""},
    {"GenSpiral3D_A", GenSpiral3D_A, METH_VARARGS, ""},
    {"GenSpiral3D_B", GenSpiral3D_B, METH_VARARGS, ""},
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