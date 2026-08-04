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
```

# Kinematics of continuum body

The motion of a continuum body is a **continuous** time sequence of displacements. Thus, the material body will occupy **different configurations** at different times so that a particle occupies a series of points in space which describe a **pathline**. There is **continuity** during deformation or motion of a continuum body in the sense that:

- The material points forming a closed curve at any instant will always form a closed curve at any subsequent time.
- The material points forming a closed surface at any instant will always form a closed surface at any subsequent time and the matter within the closed surface will always remain within

```{image} figures/displacement_intro.png
:alt: Displacement
:class: bg-primary mb-1
:width: 75%
:align: center
```

## Kinematics: deformation and motion

It is convenient to identify a **reference configuration or initial condition** which all subsequent **deformed configurations** are
referenced from. Often, the configuration at $t=0$ is considered the reference configuration.

```{image} figures/displacement_simple.png
:alt: Displacement simplified
:class: bg-primary mb-1
:width: 50%
:align: center
```

The components $x_i$ of the position vector $\mathbf{x}$ of a particle, taken with respect to the reference
configuration, are called the **material or reference coordinates**.

```{image} figures/displacement_simple.png
:alt: Displacement simplified
:class: bg-primary mb-1
:width: 50%
:align: center
```

The displacement of first point is decribed as:

$$\mathbf{x} = \mathbf{X} + \mathbf{u}(\mathbf{X})$$

while displacement of second surrounding point is described as:

$$\mathbf{x} + d\mathbf{x}  = \mathbf{X} + d\mathbf{X} + \mathbf{u}(\mathbf{X} + d\mathbf{X})$$

Substituting first equation into second we got:

$$\mathbf{X} + \mathbf{u}(\mathbf{X}) + d\mathbf{x}  = \mathbf{X} + d\mathbf{X} + \mathbf{u}(\mathbf{X} + d\mathbf{X})$$

which simplifies to:

$$ d\mathbf{x} = d\mathbf{X} + \mathbf{u}(\mathbf{X} + d\mathbf{X}) - \mathbf{u}(\mathbf{X})$$

```{admonition} Taylor's theorem
:class: tip
Taylor's theorem states that any function that is infinitely differentiable may be represented by a Taylor series expansion:

$$f(X+dX)=f(X)+{\frac  {f^{\prime }(X)}{1!}}dX+{\frac  {f^{{\prime \prime }}(X)}{2!}}dX^{2}+...=\sum _{{k=0}}^{{\infty }}{\frac  {f^{{(k)}}(X)}{k!}}dX^{{k}}$$

than 

$$f(X+dX) - f(X) = {\frac  {f^{\prime }(X)}{1!}}dX+{\frac  {f^{{\prime \prime }}(X)}{2!}}dX^{2}+...=\sum _{{k=0}}^{{\infty }}{\frac  {f^{{(k)}}(X)}{k!}}dX^{{k}}$$

neglecting higher terms as $\left | d\mathbf{X} \right | \ll 1$ as $dX^{{k}}$ is very small (we explore infinitesimal volume), it is:

$$ f(X+dX) - f(X) \approx  {\frac  {f^{\prime }(X)}{1!}}dX $$
```

Similarily (for details you have to dig into your math classes notes), for vector-valued functions we can write:

$$\mathbf{u}(\mathbf{X} + d\mathbf{X}) - \mathbf{u}(\mathbf{X}) \approx  \boldsymbol{J}(\mathbf{u})d\mathbf{X}$$

where $\boldsymbol{J}(\mathbf{u})$ is *Jacobian matrix* and in strain analysis, we usually called **displacement gradient** and we use symbol $\boldsymbol{\nabla u}$. Using that for infinitesimal deformation equation:

$$d\mathbf{x} = d\mathbf{X} + \mathbf{u}(\mathbf{X} + d\mathbf{X}) - \mathbf{u}(\mathbf{X})$$

it could be written in terms of gradient as:

$$d\mathbf{x} = d\mathbf{X} + (\boldsymbol{\nabla u})d\mathbf{X}$$

where $\boldsymbol{\nabla u}$ is gradient of displacement field or **displacement gradient**.

## Displacement gradient

The **displacement gradient** is the matrix of all first-order partial derivatives of each component of the element displacement $d\mathbf{u}$ with respect to each component of the reference element $d\mathbf{X}$:

$$\boldsymbol{\nabla u} =  u_{i,j}  =  \frac{\partial u_i}{\partial X_j} =
\begin{bmatrix}
\frac{\partial u_1}{\partial X_1} & \frac{\partial u_1}{\partial X_2} & \frac{\partial u_1}{\partial X_3} \\
\frac{\partial u_2}{\partial X_1} & \frac{\partial u_2}{\partial X_2} & \frac{\partial u_2}{\partial X_3} \\
\frac{\partial u_3}{\partial X_1} & \frac{\partial u_3}{\partial X_2} & \frac{\partial u_3}{\partial X_3}
\end{bmatrix}$$

and characterise the local change of the displacement field at a material point with position vector $\mathbf{X}$. Knowing that:

$$d\mathbf{u} = d\mathbf{x} - d\mathbf{X}$$

it could be also written as:

$$d\mathbf{u} = (\boldsymbol{\nabla u})d\mathbf{X}$$


## Deformation gradient

Recalling that $d\mathbf{u} = d\mathbf{x} - d\mathbf{X}$

$$(\boldsymbol{\nabla u}) = \frac{\partial u_i}{\partial X_j} = \frac{\partial (x_i - X_i)}{\partial X_j} = \frac{\partial x_i}{\partial X_j} - \frac{\partial X_i}{\partial X_j} = \boldsymbol{F} - \boldsymbol{I}$$

where $\boldsymbol{F}$ is so called **deformation gradient**, i.e the derivative of each component of the deformed linear element $d\mathbf{x}$ with respect to each component of the reference element $d\mathbf{X}$:

$$ \boldsymbol{F} =  x_{i,j}  =  \frac{\partial x_i}{\partial X_j} =
\begin{bmatrix}
\frac{\partial x_1}{\partial X_1} & \frac{\partial x_1}{\partial X_2} & \frac{\partial x_1}{\partial X_3} \\
\frac{\partial x_2}{\partial X_1} & \frac{\partial x_2}{\partial X_2} & \frac{\partial x_2}{\partial X_3} \\
\frac{\partial x_3}{\partial X_1} & \frac{\partial x_3}{\partial X_2} & \frac{\partial x_3}{\partial X_3}
\end{bmatrix}$$

and characterizes the local deformation at a material point with position vector $\mathbf{X}$, assuming continuity. Knowing that:

$$d\mathbf{x} - d\mathbf{X} = d\mathbf{u} = (\boldsymbol{\nabla u})d\mathbf{X} = (\boldsymbol{F} - \boldsymbol{I})d\mathbf{X} = \boldsymbol{F}d\mathbf{X} - d\mathbf{X}$$

it could be also written as:

$$d\mathbf{x} = \boldsymbol{F}d\mathbf{X}$$

## Properties of deformation gradient

**Deformation gradient** $\boldsymbol{F}$ contains all the required local information about the changes in length, volumes and angles due to the deformation as follows:

- When vector $\mathbf{X}$ in the reference configuration is deformed into the vector $\mathbf{x}$, these vectors are related as: $\mathbf{x} = \boldsymbol{F}\mathbf{X}$
- The Jacobian of the deformation gradient is equal to the ratio between the local volume of the deformed configuration to the local volume in the reference configuration i.e. volume change: $J = \frac{dV}{dV_0} = \det({\boldsymbol{F})}$
- Two infinitesimal areas with $da$ and $dA$ being their magnitudes and $\mathbf{n}$ and $\mathbf{N}$ are unit vectors perpendicular to them, then the relationship is given by: $(da)\mathbf{n} = \det(\boldsymbol{F})(dA)\boldsymbol{F}^{-T}\mathbf{N}$
- An isochoric deformation is a deformation preserving local volume, i.e., $\det({\boldsymbol{F}})=1$
- A deformation is called homogeneous if $\boldsymbol{F}$ is constant at every point. Otherwise, the deformation is called non-homogeneous
- The physical restriction of possible deformation: $\det({\boldsymbol{F}}) > 0$

## Infinitesimal strain vs. finite strain

Everything above — $\boldsymbol{F}$, $\boldsymbol{\nabla u}$ — is exact, with no restriction on how large the deformation is. The later chapters build **finite strain** measures directly from $\boldsymbol{F}$ (the Green-Lagrange strain $\boldsymbol{E}=\frac{1}{2}(\boldsymbol{C}-\boldsymbol{I})=\frac{1}{2}(\boldsymbol{F}^T\boldsymbol{F}-\boldsymbol{I})$, and the Cauchy-Green tensors $\boldsymbol{B}$, $\boldsymbol{C}$), which remain valid however large the rotation or stretch is.

There is, however, a much older and simpler measure that you will meet in classical elasticity, seismology, and geodesy texts: the **infinitesimal (small) strain tensor**. It comes from splitting the displacement gradient into its symmetric and antisymmetric parts, which is always possible for any square matrix:

$$\boldsymbol{\nabla u} = \underbrace{\tfrac{1}{2}\left(\boldsymbol{\nabla u}+(\boldsymbol{\nabla u})^T\right)}_{\boldsymbol{\varepsilon}\ \text{(infinitesimal strain)}} + \underbrace{\tfrac{1}{2}\left(\boldsymbol{\nabla u}-(\boldsymbol{\nabla u})^T\right)}_{\boldsymbol{\omega}\ \text{(infinitesimal rotation)}}$$

$\boldsymbol{\varepsilon}$ is symmetric (pure shape/size change), $\boldsymbol{\omega}$ is antisymmetric (a small rigid rotation, no shape change). Expanding the Green-Lagrange strain in terms of $\boldsymbol{\nabla u}$ shows exactly how $\boldsymbol{\varepsilon}$ relates to it:

$$\boldsymbol{E} = \tfrac{1}{2}\left(\boldsymbol{F}^T\boldsymbol{F}-\boldsymbol{I}\right) = \tfrac{1}{2}\left((\boldsymbol{I}+\boldsymbol{\nabla u})^T(\boldsymbol{I}+\boldsymbol{\nabla u})-\boldsymbol{I}\right) = \boldsymbol{\varepsilon} + \tfrac{1}{2}(\boldsymbol{\nabla u})^T\boldsymbol{\nabla u}$$

This identity holds exactly, for *any* deformation. The infinitesimal strain tensor $\boldsymbol{\varepsilon}$ is a good approximation to the finite strain $\boldsymbol{E}$ only when the quadratic term $\tfrac{1}{2}(\boldsymbol{\nabla u})^T\boldsymbol{\nabla u}$ is negligible next to the linear term — i.e. when $|\boldsymbol{\nabla u}|\ll 1$. We can check both regimes directly:

```{code-cell} ipython3
def eps_and_E(du):
    du = np.asarray(du, dtype=float)
    I = np.eye(2)
    F = I + du
    eps = 0.5 * (du + du.T)
    E = 0.5 * (F.T @ F - I)
    return eps, E

# a small deformation: |grad u| ~ 0.01
eps, E = eps_and_E([[0.01, 0.005], [0.005, -0.005]])
print('eps:\n', eps)
print('E:\n', E)
print('max |E - eps| =', np.max(np.abs(E - eps)))
```

```{code-cell} ipython3
# the same F used in the polar-decomposition example (next chapter): grad(u) = F - I, |grad u| ~ O(1)
eps, E = eps_and_E([[0.732, -0.25], [1, -0.567]])
print('eps:\n', eps)
print('E:\n', E)
print('max |E - eps| =', np.max(np.abs(E - eps)))
```

For the small deformation, $\boldsymbol{\varepsilon}$ and $\boldsymbol{E}$ agree to within $10^{-4}$ — the linearization is excellent. For the $O(1)$ deformation, they disagree by up to $0.77$, comparable to the strain itself — the linearization is useless. **This is the whole reason structural geology needs finite-strain theory** (this chapter onward): natural deformations routinely involve stretches and rotations far too large for the infinitesimal approximation, so we work directly with $\boldsymbol{F}$, $\boldsymbol{B}$, $\boldsymbol{C}$ rather than their small-strain linearization $\boldsymbol{\varepsilon}$.
