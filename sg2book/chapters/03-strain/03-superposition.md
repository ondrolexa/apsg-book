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

```{code-cell} ipython3
:tags: [remove-input]
import numpy as np
import matplotlib.pyplot as plt
from apsg import *
from strain2d import plot_defgrad
```

# Superposition of deformation

$\boldsymbol{F}$ maps any undeformed vector into its deformed state. This vector can also be a position vector of a point. Therefore $\boldsymbol{F}$ also maps any point into its new position after deformation. Considering two successive deformations $\boldsymbol{F_1}$ and $\boldsymbol{F_2}$ write transformation equation....

$$\mathbf{x}_1 = \boldsymbol{F_1} \cdot \mathbf{X}$$

$$\mathbf{x}_2 = \boldsymbol{F_2} \cdot \mathbf{x}_1$$

Substitute first equation to second gives:

$$\mathbf{x}_2 = \boldsymbol{F_2} \cdot \boldsymbol{F_1} \cdot \mathbf{X}$$

so

$$\mathbf{x}_2 = \boldsymbol{F} \cdot \mathbf{X}$$

where

$$\boldsymbol{F} = \boldsymbol{F_2} \cdot \boldsymbol{F_1}$$

The total deformation gradient can be written as the product of two partial deformation gradients, where order from right to left corresponds to superposition of deformations.

```{code-cell} ipython3
F1 = defgrad2([[1, 1], [0, 1]])
F2 = defgrad2([[2 ** 0.5, 0], [0, 0.5 ** 0.5]])
plot_defgrad(F2 @ F1, title='F2⋅F1')
```

```{code-cell} ipython3
plot_defgrad(F1 @ F2, title='F1⋅F2')
```
