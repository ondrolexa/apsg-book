---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Foreword

The power of linear algebra (also referred to as matrix algebra) as an
analytical tool in science and engineering has long been accepted. The
methods developed have been exploited with increasing frequency in
recent years because of the ease of access to desktop computers.
Geologists, traditionally slow in applying numerical methods to their
subject have nonetheless made much progress in the application of these
methods. Applications have covered the fields of structural geology;
petrology (in particular metamorphic petrology); petroleum geology
(well-log analysis); computer modelling and geostatistics.

In structural geology the use of transformation matrices has found
application in the analysis of strain, while in petrology the linear nature
of chemical equations has been exploited allowing data to be manipulated
according to systems of end-members, whether they be mineral species or
oxides of their component elements. Elsewhere, as for example in studies
of continental drift, transformations of co-ordinate systems based on
traditional algebraic and trigonometrical methods have been rewritten in
matrix form, simplifying calculations as well as aiding understanding.
Linear algebra has proved to be a particularly powerful tool in multivariate
statistics, which has led to its application in analysing multicomponent
geochemical data, which is so readily available from modern analytical
instruments.

All codes in book use following imports:

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
```

```{tableofcontents}
```
