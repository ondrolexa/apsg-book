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

The components $x_i$ of the position vector $\vec{x}$ of a particle, taken with respect to the reference
configuration, are called the **material or reference coordinates**.

```{image} figures/displacement_simple.png
:alt: Displacement simplified
:class: bg-primary mb-1
:width: 50%
:align: center
```

The displacement of first point is decribed as:

$$\vec{x} = \vec{X} + \vec{u}(\vec{X})$$

while displacement of second surrounding point is described as:

$$\vec{x} + d\vec{x}  = \vec{X} + d\vec{X} + \vec{u}(\vec{X} + d\vec{X})$$

Substituting first equation into second we got:

$$\vec{X} + \vec{u}(\vec{X}) + d\vec{x}  = \vec{X} + d\vec{X} + \vec{u}(\vec{X} + d\vec{X})$$

which simplifies to:

$$ d\vec{x} = d\vec{X} + \vec{u}(\vec{X} + d\vec{X}) - \vec{u}(\vec{X})$$

```{admonition} Taylor's theorem
:class: tip
Taylor's theorem states that any function that is infinitely differentiable may be represented by a Taylor series expansion:

$$f(X+dX)=f(X)+{\frac  {f^{\prime }(X)}{1!}}dX+{\frac  {f^{{\prime \prime }}(X)}{2!}}dX^{2}+...=\sum _{{k=0}}^{{\infty }}{\frac  {f^{{(k)}}(X)}{k!}}dX^{{k}}$$

than 

$$f(X+dX) - f(X) = {\frac  {f^{\prime }(X)}{1!}}dX+{\frac  {f^{{\prime \prime }}(X)}{2!}}dX^{2}+...=\sum _{{k=0}}^{{\infty }}{\frac  {f^{{(k)}}(X)}{k!}}dX^{{k}}$$

neglecting higher terms as $\left | d\vec{X} \right | \ll 1$ as $dX^{{k}}$ is very small (we explore infinitesimal volume), it is:

$$ f(X+dX) - f(X) \approx  {\frac  {f^{\prime }(X)}{1!}}dX $$
```

Similarily (for details you have to dig into your math classes notes), for vector-valued functions we can write:

$$\vec{u}(\vec{X} + d\vec{X}) - \vec{u}(\vec{X}) \approx  \boldsymbol{J}(\vec{u})d\vec{X}$$

where $\boldsymbol{J}(\vec{u})$ is *Jacobian matrix* and in strain analysis, we usually called **displacement gradient** and we use symbol $\boldsymbol{\nabla u}$. Using that for infinitesimal deformation equation:

$$d\vec{x} = d\vec{X} + \vec{u}(\vec{X} + d\vec{X}) - \vec{u}(\vec{X})$$

it could be written in terms of gradient as:

$$d\vec{x} = d\vec{X} + (\boldsymbol{\nabla u})d\vec{X}$$

where $\boldsymbol{\nabla u}$ is gradient of displacement field or **displacement gradient**.

## Displacement gradient

The **displacement gradient** is the matrix of all first-order partial derivatives of each component of the element displacement $d\vec{u}$ with respect to each component of the reference element $d\vec{X}$:

$$\boldsymbol{\nabla u} =  u_{i,j}  =  \frac{\partial u_i}{\partial X_j} =
\begin{bmatrix}
\frac{\partial u_1}{\partial X_1} & \frac{\partial u_1}{\partial X_2} & \frac{\partial u_1}{\partial X_3} \\
\frac{\partial u_2}{\partial X_1} & \frac{\partial u_2}{\partial X_2} & \frac{\partial u_2}{\partial X_3} \\
\frac{\partial u_3}{\partial X_1} & \frac{\partial u_3}{\partial X_2} & \frac{\partial u_3}{\partial X_3}
\end{bmatrix}$$

and characterise the local change of the displacement field at a material point with position vector $\vec{X}$. Knowing that:

$$d\vec{u} = d\vec{x} - d\vec{X}$$

it could be also written as:

$$d\vec{u} = (\boldsymbol{\nabla u})d\vec{X}$$


## Deformation gradient

Recalling that $d\vec{u} = d\vec{x} - d\vec{X}$

$$(\boldsymbol{\nabla u}) = \frac{\partial u_i}{\partial X_j} = \frac{\partial (x_i - X_i)}{\partial X_j} = \frac{\partial x_i}{\partial X_j} - \frac{\partial X_i}{\partial X_j} = \boldsymbol{F} - \boldsymbol{I}$$

where $\boldsymbol{F}$ is so called **deformation gradient**, i.e the derivative of each component of the deformed linear element $d\vec{x}$ with respect to each component of the reference element $d\vec{X}$:

$$ \boldsymbol{F} =  x_{i,j}  =  \frac{\partial x_i}{\partial X_j} =
\begin{bmatrix}
\frac{\partial x_1}{\partial X_1} & \frac{\partial x_1}{\partial X_2} & \frac{\partial x_1}{\partial X_3} \\
\frac{\partial x_2}{\partial X_1} & \frac{\partial x_2}{\partial X_2} & \frac{\partial x_2}{\partial X_3} \\
\frac{\partial x_3}{\partial X_1} & \frac{\partial x_3}{\partial X_2} & \frac{\partial x_3}{\partial X_3}
\end{bmatrix}$$

and characterizes the local deformation at a material point with position vector $\vec{X}$, assuming continuity. Knowing that:

$$d\vec{x} - d\vec{X} = d\vec{u} = (\boldsymbol{\nabla u})d\vec{X} = (\boldsymbol{F} - \boldsymbol{I})d\vec{X} = \boldsymbol{F}d\vec{X} - d\vec{X}$$

it could be also written as:

$$d\vec{x} = \boldsymbol{F}d\vec{X}$$

## Properties of deformation gradient

**Deformation gradient** $\boldsymbol{F}$ contains all the required local information about the changes in length, volumes and angles due to the deformation as follows:

- When vector $\vec{X}$ in the reference configuration is deformed into the vector $\vec{x}$, these vectors are related as: $\vec{x} = \boldsymbol{F}\vec{X}$
- The Jacobian of the deformation gradient is equal to the ratio between the local volume of the deformed configuration to the local volume in the reference configuration i.e. volume change: $J = \frac{dV}{dV_0} = \det({\boldsymbol{F})}$
- Two infinitesimal areas with $da$ and $dA$ being their magnitudes and $\vec{n}$ and $\vec{N}$ are unit vectors perpendicular to them, then the relationship is given by: $(da)\vec{n} = \det(\boldsymbol{F})(dA)\boldsymbol{F}^{-T}\vec{N}$
- An isochoric deformation is a deformation preserving local volume, i.e., $\det({\boldsymbol{F}})=1$
- A deformation is called homogeneous if $\boldsymbol{F}$ is constant at every point. Otherwise, the deformation is called non-homogeneous
- The physical restriction of possible deformation: $\det({\boldsymbol{F}}) > 0$
