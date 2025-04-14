# RETM (Reflectance and Emittance Thermophysical Model) or TREM (Thermophysical Reflectance and Emittance Model)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developer
- Jin BENIYAMA [mail](mailto:jinbeniyama@gmail.com)
- hoge
- hoge

## Overview
A model that implements the theory of reflectance and emittance spectroscopy.
Under heavy development by J.B. etc.

## Structure
```
retm_or_trem/
  README.md
  docs/
  notebooks/
  retm/ 
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
