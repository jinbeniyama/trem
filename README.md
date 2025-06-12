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

## Test (in /tests)
### Perform TPM
A test TPM can be performed with a following command.
We need to understand what `obj file`, `spin file`, `observation file`, and `ephemeris file` are.
Please refer to `readme.pdf` (the document of Marco's code), [Documentation of DAMIT](https://astro.troja.mff.cuni.cz/projects/damit/), `src/make_obseph.py` for the details of each file.
```
# Do TPM
bash TPM_test.sh SHAPE.obj spin.txt WISE.obs eph.txt
```

### Make observation and ephemeris files of Ryugu
The observation file and ephemeris file can be made using `src/make_obseph.py`.
We need to prepare observational results as `flux_test.txt` beforehand.
The light-time correction is not performed by default, but can be performed with `--ltcor` option.

```
# Make observation and ephemeris files of Ryugu
make_obseph.py Ryugu flux_test.txt --out_obs Ryugu_obs.txt --out_eph Ryugu_eph.txt
# Make observation and ephemeris files of Ryugu w/light-time correction
make_obseph.py Ryugu flux_test.txt --out_obs Ryugu_obs.txt --out_eph Ryugu_eph.txt --ltcor
```

## Dependencies
This repository is depending on `Python`, `NumPy`, `pandas`, `SciPy`, `Astropy`, `Astroquery`.

## LICENCE
This software is released under the MIT License.
