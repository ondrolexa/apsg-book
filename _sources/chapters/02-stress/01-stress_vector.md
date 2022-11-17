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

# Stress vector

Consider two blocks of different cross-sections. Intuitively, the blocks whose cross-section is smaller is going to deform a lot
more than the other. While in rigid body mechanics, the concept of force is sufficient to describe or predict the motion of the body, in deformable bodies it is not. 

```{image} figures/twiss_stress.png
:alt: Stress
:class: bg-primary mb-1
:width: 50%
:align: center
```

We will need to understand the concept of **stress**. The term **stress** ($\sigma$) is used to express the loading in terms of force applied to a certain cross-sectional area of an object.

From the perspective of loading, stress is the applied force or system of forces that tends to deform a body.

From the perspective of what is happening within a material, stress is the internal distribution of forces within a body that balance and react to the loads applied to it. The stress distribution may or may not be uniform, depending on the nature of the loading condition. 

```{image} figures/simple_stress.png
:alt: Simple stress
:class: bg-primary mb-1
:width: 50%
:align: center
```

## Traction

Simplifying assumptions are often used to represent force acting on area as **traction vector** - simply the force vector divided by that area. 

$$\boldsymbol{T} = \frac{\boldsymbol{F}}{A}$$

$\boldsymbol{T}$ has units of stress, but as it is a vector, all the usual rules for vectors apply to it.

The unit of stress is Pascal [Pa] or Newtons per square meter [$\frac{N}{m^2}$]. A stress of 1 Pa is very small, for example, the load due to 1m of water is about 10$^4$ Pa. So, in geology we often use metric prefixes:
 - 1 hPa = 10$^2$ Pa (hektopascal)
 - 1 kPa = 10$^3$ Pa (kilopascal)
 - 1 MPa = 10$^6$ Pa (megapascal)
 - 1 GPa = 10$^9$ Pa (gigapascal)

Other units are bars and atm: 1 MPa = 10$^6$ Pa = 10 bars = 9.869232667 atm

The **traction** is a **vector quantity** that acts at a point on an imaginary or real surface of arbitrary orientation.

```{admonition} Cauchy reciprocal theorem
:class: tip
The traction at a point on a surface is equal and opposite to the traction that at that same point for the same surface with opposite outward unit normal vector.
```

```{image} figures/tractions.png
:alt: Traction
:class: bg-primary mb-1
:width: 50%
:align: center
```

## Stress components

In general, a stress acting on a plane represented by traction $\boldsymbol{T}^{(\boldsymbol{n})}$ may be expressed as a sum
of shear and normal components.

```{image} figures/stress_onplane.png
:alt: Stress on plane
:class: bg-primary mb-1
:width: 60%
:align: center
```

```{note}
**Normal stress**: The component of stress acting perpendicular to the plane

**Shear stress**: The component of stress acting parallel to the plane
```

The magnitude of the normal stress component $\sigma _n$ of any stress vector $\boldsymbol{T}^{(\boldsymbol{n})}$ acting on an arbitrary
plane with normal unit vector $n$ at a given point, in terms of the stress tensor $\boldsymbol{\sigma}$, could be calculated a scalar projection of the stress vector onto the normal unit vector:

$$\left\| {{{\vec \sigma }_n}} \right\| = {{\bf{T}}^{\left( {\bf{n}} \right)}} \cdot \bf{n}$$

$${{\vec \sigma }_n} = \left\| {{{\vec \sigma }_n}} \right\|\bf{n}$$

The magnitude of the shear stress component $\tau_n$, acting in the plane spanned by the two vectors $\boldsymbol{T}^{(\boldsymbol{n})}$ and $n$, can then be found using the Pythagorean theorem:

$$\left\| {{{\vec \tau }_n}} \right\| = \sqrt {{{\left\| {{{\bf{T}}^{\left( {\bf{n}} \right)}}} \right\|}^2} - {{\left\| {{{\vec \sigma }_n}} \right\|}^2}}$$

$${{\vec \tau }_n} = {{\bf{T}}^{\left( {\bf{n}} \right)}} - {{\vec \sigma }_n}$$


## Stress at point
It must be noted that the stresses in most 2-D or 3-D solids are actually more complex and need be defined more methodically.

The internal force acting on a small area of a plane represented by traction vector can be resolved into normal and shear stress components. These stresses are average stresses as the area is finite, but when the area is allowed to approach zero, the stresses become stresses at a point.

Since stresses are defined in relation to the plane that passes through the point under consideration, and the number of such planes is infinite, there appear an infinite set of stresses at a point.

Fortunately, it can be proven that the stresses on any plane can be computed from the stresses on three orthogonal planes passing through the point... and **tensors** will help us.
