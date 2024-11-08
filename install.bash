#!/bin/bash

pip install numpy
pip install matplotlib

python -m build .
pip install .
