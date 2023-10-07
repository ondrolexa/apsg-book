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
from apsg import *
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

$$\mathbf{t} = \frac{\mathbf{f}}{A}$$

$\mathbf{t}$ has units of stress, but as it is a vector, all the usual rules for vectors apply to it.

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

In general, a stress acting on a plane represented by traction $\mathbf{t}$ may be expressed as a sum
of shear and normal components.

$$\mathbf{t} = \mathbf{\sigma}_n + \mathbf{\tau}_{x_1} + \mathbf{\tau}_{x_1} = \mathbf{\sigma}_n + \mathbf{\tau}_n$$

```{image} figures/stress_onplane.png
:alt: Stress on plane
:class: bg-primary mb-1
:width: 60%
:align: center
```

```{note}
**Normal stress**: The component of stress acting perpendicular to the plane ($\mathbf{\sigma}_n$)

**Shear stress**: The component of stress acting parallel to the plane ($\mathbf{\tau}_n$)
```

The magnitude of the normal stress component $\mathbf{\sigma}_n$ of any stress vector $\mathbf{t}$ acting on an arbitrary
plane with normal unit vector $\mathbf{n}$ at a given point, could be calculated a scalar projection of the stress vector
onto the normal unit vector:

$$\left\| \mathbf{\sigma}_n \right\| = \mathbf{t} \cdot \mathbf{n}$$

$$\mathbf{\sigma}_n = \left\| \mathbf{\sigma}_n \right\|\mathbf{n}$$

The magnitude of the shear stress component $\mathbf{\tau}_n$, acting in the plane spanned by the two vectors $\mathbf{t}$ and
$\mathbf{\sigma}_n$, can then be found using the Pythagorean theorem:

$$\left\| \mathbf{\tau} \right\| = \sqrt {{\left\| \mathbf{t} \right\|}^2 - {\left\| \mathbf{\sigma}_n \right\|}^2}$$

$$\mathbf{\tau}_n = \mathbf{t} - \mathbf{\sigma}_n$$

```{code-cell} ipython3
n = fol(150, 60)
T = 5 * vec('x')  # Traction vector oriented parallel to x-axis with magnitude 5
T.dot(n)  # normal stress magnitude
```

```{code-cell} ipython3
sn = T.proj(n)  # project tranction vector T onto normal n - normal stress
abs(sn)  # magnitude of normal stress
```

```{code-cell} ipython3
tau = T - sn  # shear stress
abs(tau)  # magnitude of shear stress
```

or we can use vector rejection:

```{code-cell} ipython3
tau = T.reject(n)  # shear stress
abs(tau)  # magnitude of shear stress
```

## Stress at point
It must be noted that the stresses in most 2-D or 3-D solids are actually more complex and need be defined more methodically.

The internal force acting on a small area of a plane represented by traction vector can be resolved into normal and shear stress components. These stresses are average stresses as the area is finite, but when the area is allowed to approach zero, the stresses become stresses at a point.

Since stresses are defined in relation to the plane that passes through the point under consideration, and the number of such planes is infinite, there appear an infinite set of stresses at a point.

Fortunately, it can be proven that the stresses on any plane can be computed from the stresses on three orthogonal planes passing through the point. This will bring us to the definition of the **stress tensor**...
