from setuptools import setup, Extension
import numpy

module = Extension(
    'mrtrjgen.ext', 
    sources = ['./mrtrjgen/ext/main.cpp', './mrtrjgen/ext/Spiral3D.cpp'],
    include_dirs = ["mrtrjgen/ext/", numpy.get_include()],
    language = 'c++'
    )

setup(
    name = 'mrtrjgen', 
    py_modules = ["mrtrjgen"],
    ext_modules = [module],
    packages=["mrtrjgen"],
    package_data={"mrtrjgen" : ["./ext/*.h"]}
    )
