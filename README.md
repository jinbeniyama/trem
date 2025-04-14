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
retm_or_trem/
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
This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.


## LICENCE
This software is released under the MIT License.
