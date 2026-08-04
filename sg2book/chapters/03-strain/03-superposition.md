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

The two resulting ellipses are not the same — applying $\boldsymbol{F}_1$ then $\boldsymbol{F}_2$ gives a different final shape and orientation than applying $\boldsymbol{F}_2$ then $\boldsymbol{F}_1$, because matrix multiplication does not commute in general: $\boldsymbol{F}_2\cdot\boldsymbol{F}_1 \ne \boldsymbol{F}_1\cdot\boldsymbol{F}_2$.

## Coaxial and non-coaxial deformation

This non-commutativity is the mathematical signature of a fundamental distinction in structural geology:

```{admonition} Coaxial vs. non-coaxial deformation
:class: tip
A progressive deformation is **coaxial** if the principal strain axes of every increment stay parallel to the principal axes of the finite (total) strain throughout the deformation history — successive increments accumulate without rotating relative to each other. **Pure shear** is the classic coaxial deformation.

A progressive deformation is **non-coaxial** if the principal axes of successive increments rotate relative to the finite strain axes as the deformation proceeds — each increment stretches a different orientation than the one before it. **Simple shear** is the classic non-coaxial deformation.

Two deformation increments that share the same principal directions (i.e. are simultaneously diagonal in some common basis) always commute; increments that do not share principal directions generally do not. Non-commutativity of $\boldsymbol{F}_1$ and $\boldsymbol{F}_2$ is therefore direct evidence that the corresponding progressive deformation is non-coaxial.
```

We can check the coaxial case directly: two stretches sharing the same (here, coordinate-aligned) principal axes commute, unlike the simple-shear/pure-shear pair above:

```{code-cell} ipython3
G1 = defgrad2([[2, 0], [0, 0.5]])
G2 = defgrad2([[1.5, 0], [0, 1 / 1.5]])
G1 @ G2 == G2 @ G1  # both diagonal in the same basis -> coaxial -> commute
```

We formalize pure shear and simple shear in terms of the polar decomposition in the next chapter.
