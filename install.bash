#!/bin/bash

pip install numpy
pip install matplotlib
pip install build

python -m build .
pip install .
