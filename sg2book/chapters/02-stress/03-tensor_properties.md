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
```

# Stress tensor properties

## Transformation rule of the stress tensor

From an $x_i$ - system to an $x_i'$ - system, the components $\sigma_{ij}$ in the initial system are transformed into the components $\sigma_{ij}'$ in the new system according to the **tensor transformation rule**:

$$\boldsymbol{{\sigma}'} = \boldsymbol{R} \boldsymbol{\sigma} \boldsymbol{R}^{T}$$

Note that inverse rule is: $\boldsymbol{{\sigma}} = \boldsymbol{R}^{T} \boldsymbol{{\sigma}'} \boldsymbol{R}$

```{figure} figures/Stress_transformation_3D.png
:name: fig-stress-transform-3d
:alt: Stress transformation in 3D
:class: bg-primary mb-1
:width: 75%
:align: center

Stress tensor transformation between coordinate systems related by rotation $\boldsymbol{R}$.
```

where $\boldsymbol{R}$ is a rotation matrix with components $a_{ij}$ ({numref}`fig-stress-transform-3d`). In matrix form this is:

$$\left[{{\begin{matrix}\sigma '_{{11}}&\sigma '_{{21}}&\sigma '_{{31}}\\\sigma '_{{12}}&\sigma '_{{22}}&\sigma '_{{32}}\\\sigma '_{{13}}&\sigma '_{{23}}&\sigma '_{{33}}\\\end{matrix}}}\right]=\left[{{\begin{matrix}a_{{11}}&a_{{12}}&a_{{13}}\\a_{{21}}&a_{{22}}&a_{{23}}\\a_{{31}}&a_{{32}}&a_{{33}}\\\end{matrix}}}\right]\left[{{\begin{matrix}\sigma _{{11}}&\sigma _{{21}}&\sigma _{{31}}\\\sigma _{{12}}&\sigma _{{22}}&\sigma _{{32}}\\\sigma _{{13}}&\sigma _{{23}}&\sigma _{{33}}\\\end{matrix}}}\right]\left[{{\begin{matrix}a_{{11}}&a_{{21}}&a_{{31}}\\a_{{12}}&a_{{22}}&a_{{32}}\\a_{{13}}&a_{{23}}&a_{{33}}\\\end{matrix}}}\right]$$

In following example we will define stress tensor in orientation given by foliation and lineation plane:

```{code-cell} ipython3
S = stress.from_comp(xx=8, yy=6, zz=2)  # positive = compression (see sign convention, previous chapter)
p = pair(90, 60, 90, 60)
R = rotation.from_pair(p)  # rotation to local coordinates
lin(0, 0).transform(R)
```

```{code-cell} ipython3
fol(0, 0).transform(R)
```

```{code-cell} ipython3
Sr = S.transform(R)
abs(Sr.cauchy(lin(90, 60)))  # magnitude of stress on plane perpendicular to 90/60
```

```{code-cell} ipython3
Sr = S.transform(R)
abs(Sr.cauchy(fol(90, 60)))  # magnitude of stress on plane 90/60
```

## Principal stresses and principal directions

At every point in a stressed body there are at least three planes, called **principal planes**, with normal vectors $\boldsymbol{n}$, called **principal directions**, where the corresponding stress vector is perpendicular to the plane, i.e., parallel or in the same
direction as the normal vector $\boldsymbol{n}$, and where there are no normal shear stresses $\tau_{\mathrm{n}}$. The three stresses normal to these principal planes are called **principal stresses**.

Note three orientations with zero shear stress in following example. Their positions are **principal directions**:

```{code-cell} ipython3
p = pair(150, 60, 90, 41)
S2 = S.transform(rotation.from_pair(p))
s = StereoNet()
s.grid.apply_func(S2.shear_stress)
s.contour(levels=10)
s.show()
```

```{code-cell} ipython3
S2.eigenlins()  # principal directions
```

```{code-cell} ipython3
S2.eigenfols()  # principal planes
```

The components $\sigma_{ij}$ of the stress tensor depend on the orientation of the coordinate system at the point under consideration. However, the stress tensor itself is a physical quantity and as such, it is independent of the coordinate system chosen to represent it. There are certain **invariants** associated with every tensor which are also independent of the coordinate system. One set of such invariants are the **principal stresses** of the stress tensor, which are just the **eigenvalues** of the stress tensor. Their direction vectors are the **principal directions** or **eigenvectors**.

```{code-cell} ipython3
S.eigenvalues()  # principal stresses of S
```

```{code-cell} ipython3
S2.eigenvalues()  # principal stresses of S2
```

## Mohr diagram in 3D

The single 2D Mohr circle from the previous chapter generalizes to **three** circles in 3D, one for each pair of principal stresses. For principal stresses $\sigma_1 \ge \sigma_2 \ge \sigma_3$, the $(\sigma_n,\tau)$ point for *any* plane orientation always falls **inside** the outer circle (built from $\sigma_1,\sigma_3$) and **outside** the two inner circles (built from $\sigma_1,\sigma_2$ and $\sigma_2,\sigma_3$) — the shaded region between them is exactly the set of physically realizable $(\sigma_n,\tau)$ pairs for this stress state.

We verify this the same way as before: instead of trusting the claim, sample many random plane orientations, compute their true $(\sigma_n,\tau)$ from `S2.cauchy(n)`, and check every single point falls inside the admissible band bounded by the three circles built from `S2.eigenvalues()`.

```{code-cell} ipython3
s1, s2, s3 = sorted(S2.eigenvalues(), reverse=True)

sn_vals, tau_vals = [], []
for _ in range(300):
    n = vec.random()
    t = S2.cauchy(n)
    sn_vals.append(abs(t.project(n)))
    tau_vals.append(abs(t.reject(n)))

fig, ax = plt.subplots()
for (a, b), color, lbl in [((s1, s3), 'C1', r'$\sigma_1,\sigma_3$ (outer)'),
                            ((s1, s2), 'C2', r'$\sigma_1,\sigma_2$'),
                            ((s2, s3), 'C2', r'$\sigma_2,\sigma_3$')]:
    c, r = (a + b) / 2, abs(a - b) / 2
    ax.add_patch(plt.Circle((c, 0), r, fill=False, color=color, label=lbl))
ax.plot(sn_vals, tau_vals, '.', ms=3, color='C0', alpha=0.6, label='random planes (S2.cauchy(n))')
ax.set_xlim(s3 - 1, s1 + 1)
ax.set_ylim(-0.5, (s1 - s3) / 2 + 1)
ax.set_xlabel(r'$\sigma_n$')
ax.set_ylabel(r'$\tau$')
ax.set_aspect(1)
ax.legend(fontsize=8)
ax.set_title('3D Mohr diagram for S2');
```

Every sampled point lands inside the outer circle and outside both inner circles, exactly as predicted — none of the stereonet contour values plotted above (`S2.shear_stress`) can exceed the radius of the outer circle, $(\sigma_1-\sigma_3)/2$, which is also the **maximum possible shear stress** at this point, realized on planes bisecting $\sigma_1$ and $\sigma_3$.

## Stress invariants

A stress vector parallel to the normal unit vector $\boldsymbol{n}$ is given by:

$$\boldsymbol{T}^{(\boldsymbol {n} )} = \boldsymbol{\sigma} \cdot \boldsymbol{n} = \lambda \boldsymbol{n}$$

where $\lambda$ is a constant of proportionality, and in this particular case corresponds to the magnitudes $\sigma_n$ of the normal stress vectors or principal stresses.

$$\begin{aligned}
  \boldsymbol{\sigma} \cdot \boldsymbol{n} &= \lambda \boldsymbol{n} \\
  \boldsymbol{\sigma} \cdot \boldsymbol{n} - \lambda \boldsymbol{n} &= \boldsymbol{0} \\
  \left(\boldsymbol{\sigma} - \lambda \delta _{ij}\right)\boldsymbol{n} &= \boldsymbol{0}\\
\end{aligned}$$

This is a homogeneous system (note zero right side) of linear equations where components $n_{i}$ of vector $\boldsymbol{n}$ are the unknowns.

To obtain a nontrivial (non-zero) solution for $n_{i}$, the determinant matrix of the coefficients must be equal to zero, i.e. the system is singular. Thus,

$$\left|\boldsymbol{\sigma} - \lambda \delta _{ij}\right| = -\lambda ^{3}+I_{1}\lambda ^{2}-I_{2}\lambda +I_{3}=0$$

where

$$\begin{aligned}I_{1}&=\sigma _{11}+\sigma _{22}+\sigma _{33}=\sigma _{kk}={\text{tr}}({\boldsymbol {\sigma }})\\I_{2}&=\sigma _{11}\sigma _{22}+\sigma _{22}\sigma _{33}+\sigma _{11}\sigma _{33}-\sigma _{12}^{2}-\sigma _{23}^{2}-\sigma _{31}^{2}\\&={\frac {1}{2}}\left(\sigma _{ii}\sigma _{jj}-\sigma _{ij}\sigma _{ji}\right)={\frac {1}{2}}\left[\left({\text{tr}}({\boldsymbol {\sigma }})\right)^{2}-{\text{tr}}({\boldsymbol {\sigma }}^{2})\right]\\I_{3}&=\det(\sigma _{ij})=\det({\boldsymbol {\sigma }})\end{aligned}$$

The $\sigma _{1}=\max \left(\lambda _{1},\lambda _{2},\lambda _{3}\right)$, $\sigma _{3}=\min \left(\lambda _{1},\lambda _{2},\lambda _{3}\right)$ and $\sigma _{2}=I_{1}-\sigma _{1}-\sigma _{3}$, are the principal stresses. The coefficients $I_{1}$, $I_{2}$ and $I_{3}$, are called the **first, second, and third stress invariants**.

```{code-cell} ipython3
S.I1, S.I2, S.I3  # stress invariants
```

```{code-cell} ipython3
S2.I1, S2.I2, S2.I3  # stress invariants
```

## Tensor decomposition

Often it is convenient to decompose the stress tensor into volumetric $\mathbf{V}$ (also known as hydrostatic) and deviatoric $\mathbf{S}$ (distortional) parts. Applications of such decompositions can be found in metal plasticity, soil mechanics, and biomechanics.

$$\mathbf{V} = \frac{\text{tr}(\boldsymbol{\sigma})}{3}\delta_{ij} \quad \text{and} \quad \mathbf{S} = \boldsymbol{\sigma} - \mathbf{V}$$

```{code-cell} ipython3
S2.hydrostatic
```

```{code-cell} ipython3
S2.deviatoric
```
