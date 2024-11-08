# Trajectory Library for Magnetic Resonance Imaging

This is the published repository of submitted abstract of ISMRM 2025: `A Novel 3D Spiral Trajectory for Comprehensive k-Space Coverage and Efficient Sampling`

1. Run `bash install.bash` to install the essential packages and the `mrtrjgen` package, you may need to install `python` and `pip` first.
1. Please cd to `Example` before you run `python Example_*.py`, DO NOT run example in the `root directory`.
1. If you want to see the proof of constant slew rate sampling of Spiral-3D or Seiffert-Spiral, please refers to `Example/Proof_Spiral3D_TypeA.py` and `Example/Proof_Spiral3D_TypeB.py`. These two Python scripts will save the derived coefficients of quadratic equation to `Example/latex.md` by default.