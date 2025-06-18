# TREM (Thermophysical Reflectance and Emittance Model)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developers
- Jin BENIYAMA [mail](mailto:jinbeniyama@gmail.com)
- Koki YUMOTO

## Overview
A model that implements the theory of reflectance and emittance spectroscopy.
Under heavy development by J.B. and K.Y.

## Structure
```
trem/
  README.md
  c/ # TPM code (in prep.)
  docs/
  notebooks/
  trem/ 
    emittance/   # To handle thermal emission (i.e., conventional TPM)
    reflectance/ # To handle NIR wavelength
    test/        # For test, whatever
  src/ # executable script (.py, .sh, etc.)
    main.py
    utils/
      helper.py
  tests/ # For test, whatever
```

## Dependencies
This repository is depending on `Python`, `NumPy`, `pandas`, `SciPy`, `Astropy`, `Astroquery`.

## LICENCE
This software is released under the MIT License.
