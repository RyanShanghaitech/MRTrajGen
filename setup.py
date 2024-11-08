from setuptools import setup, Extension
import numpy

module = Extension(
    'mrtrjgen.ext', 
    sources = 
    [
        './mrtrjgen/ext/main.cpp', 
        './mrtrjgen/ext/Spiral3D_A.cpp', 
        './mrtrjgen/ext/Spiral3D_B.cpp', 
        './mrtrjgen/ext/Trajectory.cpp'
    ],
    include_dirs = ["./ext/", numpy.get_include()],
    language = 'c++'
    )

setup(
    name = 'mrtrjgen',
    ext_modules = [module],
    packages = ["mrtrjgen"],
    package_data={"mrtrjgen" : ["./ext/*.h"]}
    )
