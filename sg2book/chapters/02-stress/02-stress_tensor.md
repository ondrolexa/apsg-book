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

# Stress tensor

The Cauchy tetrahedron derivation and stress tensor notation used here follow {cite}`allmendinger2012structural`.

## Cauchy formula in 2D

```{figure} figures/Cauchy_2D_triangle_v2.png
:name: fig-cauchy-2d
:alt: Cauchy in 2D
:class: bg-primary mb-1
:width: 80%
:align: center

Equilibrium of a 2D triangular element used to derive the Cauchy stress relation.
```

The equilibrium of forces, i.e. Euler’s first law of motion (Newton’s second law of motion) for the 2D triangle in {numref}`fig-cauchy-2d`:

$$\mathbf{t}A = \mathbf{t}_1 A_1 + \mathbf{t}_2 A_2$$

The area of the faces of the triangle perpendicular to the axes can be found by projecting $A$ into each face

$$\mathbf{t}A = \mathbf{t}_1 n_1 A + \mathbf{t}_2 n_2 A$$

than divided by $A$ gives:

$$\mathbf{t} = \mathbf{t}_1 n_1 + \mathbf{t}_2 n_2$$

This equation could be written by components

$$\begin{aligned}
t_1 &= \sigma_{11} n_1 + \sigma_{21} n_2 \\
t_2 &= \sigma_{12} n_1 + \sigma_{22} n_2 
\end{aligned}$$

or in matrix form

$$\begin{bmatrix}t_1\\t_2\end{bmatrix}=
\begin{bmatrix}
\sigma_{11} & \sigma_{21}\\ 
\sigma_{12} & \sigma_{22}
\end{bmatrix}
\cdot
\begin{bmatrix}
n_1\\ 
n_2
\end{bmatrix}
\quad \text{or} \quad
\mathbf{t} = \boldsymbol{\sigma} \cdot \mathbf{n}
\quad \text{where} \quad
\boldsymbol{\sigma} = \begin{bmatrix}\mathbf{t}_1 & \mathbf{t}_2\end{bmatrix}$$

Because the $\sigma_{12} = \sigma_{21}$ (conservation of angular momentum), the $\boldsymbol{\sigma}$ must be **symmetric**.

```{note}
Look closely at where each component sits in the matrix above: entry (row $i$, column $j$) holds $\sigma_{ji}$, not $\sigma_{ij}$ — the matrix is built from *columns* $\mathbf{t}_1,\mathbf{t}_2$ (the traction vectors), which places the printed layout one transpose away from the "row $i$ = plane, column $j$ = direction" index convention we state below in *Indexes of stress components*. This is harmless here only because $\boldsymbol{\sigma}=\boldsymbol{\sigma}^T$ for a symmetric tensor — $\Sigma$ and $\Sigma^T$ are numerically identical. It stops being harmless for a genuinely asymmetric tensor: the displacement gradient $\boldsymbol{\nabla u}$ and deformation gradient $\boldsymbol{F}$ used later in this book are *not* symmetric in general, so the same row/column placement would matter there and must not be pattern-matched blindly from this chapter.
```

```{admonition} Matrix as function
:class: tip
Note that $\boldsymbol{\sigma}$, known as **stress tensor** operates as a function, i.e. each normal unit vector $\mathbf{n}$ is transformed to traction vector $\mathbf{t}$ acting on the plane normal to $\mathbf{n}$. In this way, stress tensor describing stress on arbitrary plane, i.e describing *stress in point*.
```

```{note}
Throughout this book we follow the geosciences and rock-mechanics sign convention: **compression is positive, tension is negative**. This is the convention APSG uses for its `sigma1`/`sigma2`/`sigma3` shortcuts ($\sigma_1$ is the *most compressive* principal stress). See *Sign convention* below for the formal definition. The example below therefore uses positive diagonal components to represent a triaxially compressive stress state, as expected for rocks at depth.
```

```{code-cell} ipython3
S = stress2([[8, -2],
             [-2,  5]])
n = vec2(60)  # unit length vector oriented 60 degrees from x
S.cauchy(n)  # traction vector
```

## Cauchy formula in 3D

The tetrahedron ({numref}`fig-cauchy-3d`) is formed by slicing the infinitesimal element along an arbitrary plane with normal $\mathbf{n}$. The traction vector on this plane is denoted by $\mathbf{t}$. The traction vectors acting on the faces of the tetrahedron are denoted as $\mathbf{t}_1$, $\mathbf{t}_2$, and $\mathbf{t}_3$, and are by definition the components of the stress tensor $\boldsymbol{\sigma}$.

```{figure} figures/Cauchy_tetrahedron.png
:name: fig-cauchy-3d
:alt: Cauchy in 3D
:class: bg-primary mb-1
:width: 60%
:align: center

The Cauchy tetrahedron: slicing an infinitesimal element along an arbitrary plane with normal $\mathbf{n}$.
```

Similarly to 2D case, the equilibrium of forces, i.e. Euler’s first law of motion (Newton’s second law of motion), gives:

$$\mathbf{t} A = \mathbf{t}_1 A_1 + \mathbf{t}_2 A_2 + \mathbf{t}_3 A_3$$

The area of the faces of the tetrahedron perpendicular to the axes can be found by projecting $A$ into each face:

$$\mathbf{t} A = \mathbf{t}_1 n_1 A + \mathbf{t}_2 n_2 A + \mathbf{t}_3 n_3 A$$

than divided by $A$ gives:

$$\mathbf{t} = \mathbf{t}_1 n_1 + \mathbf{t}_2 n_2 + \mathbf{t}_3 n_3$$


```{admonition} What is tensor
:class: tip
**Tensors** are algebraic objects that describe linear relationship between vectors, scalars, or tensors. Here, any linear connection between two physical vector quantities is called a **tensor**, reflecting original use to describe the "tensions" in a material.

$$\left[{\begin{matrix} u_1 \\ u_2 \\ u_3 \end{matrix}}\right] = \left[{\begin{matrix}
    a_{11} & a_{21} & a_{31} \\
    a_{12} & a_{22} & a_{32} \\
    a_{13} & a_{23} & a_{33}
\end{matrix}}\right] \left[{\begin{matrix} v_1 \\ v_2 \\ v_3 \end{matrix}}\right]$$

$$\mathbf{u} = \mathbf{A} \cdot \mathbf{v}$$
```

## Cauchy stress tensor

In continuum mechanics, the Cauchy **stress tensor** $\boldsymbol{\sigma}$ is a second order tensor,
with nine components $\sigma_{ij}$ ({numref}`fig-stress-components`), that completely define the state of stress at a point inside a material.

According to the *principle of conservation of angular momentum*, equilibrium requires that the summation of moments with respect to an arbitrary point is zero, which leads to the conclusion that the stress tensor is symmetric, thus having only *six independent stress components*, instead of the original nine.

```{figure} figures/Components_stress_tensor_cartesian.png
:name: fig-stress-components
:alt: Stress tensor components
:class: bg-primary mb-1
:width: 60%
:align: center

The nine components of the Cauchy stress tensor on the faces of a Cartesian volume element.
```

$$\boldsymbol{\sigma} = \sigma_{ij} = \left[{\begin{matrix} \mathbf{t}_1 & \mathbf{t}_2 & \mathbf{t}_3 \end{matrix}}\right] =
\left[{\begin{matrix}
    \sigma _{11} & \sigma _{21} & \sigma _{31} \\
    \sigma _{12} & \sigma _{22} & \sigma _{32} \\
    \sigma _{13} & \sigma _{23} & \sigma _{33}
\end{matrix}}\right]$$

$$\left[{\begin{matrix}
    \sigma _{11} & \sigma _{21} & \sigma _{31} \\
    \sigma _{12} & \sigma _{22} & \sigma _{32} \\
    \sigma _{13} & \sigma _{23} & \sigma _{33}
\end{matrix}}\right] \equiv \left[{\begin{matrix}
    \sigma _{xx} & \sigma _{yx} & \sigma _{zx} \\
    \sigma _{xy} & \sigma _{yy} & \sigma _{zy} \\
    \sigma _{xz} & \sigma _{yz} & \sigma _{zz}
\end{matrix}}\right] \equiv \left[{\begin{matrix}
    \sigma _x & \tau _{xy} & \tau _{xz} \\
    \tau _{xy} & \sigma _y & \tau _{yz} \\
    \tau _{xz} & \tau _{yz} & \sigma _z \\
\end{matrix}}\right]$$
  
where $\sigma_{11}$, $\sigma_{22}$, $\sigma_{33}$ are normal stresses, and $\sigma_{12}$, $\sigma_{13}$, $\sigma_{21}$, $\sigma_{23}$, $\sigma_{31}$, $\sigma_{32}$ are shear stresses.

```{admonition} Indexes of stress components
:class: tip
The first index $i$ indicates that the stress acts on a plane normal to the $x_i$-axis, and the second index $j$ denotes the direction in which the stress acts.
```

## Cauchy's stress theorem

According to Cauchy’s fundamental theorem, also called **Cauchy's stress theorem**, merely by knowing the traction vectors on three mutually perpendicular planes, the stress vector on any other plane passing through that point can be found using transformation equation.

$$\mathbf{t} \equiv \left[{\begin{matrix} \mathbf{t}_1 & \mathbf{t}_2 & \mathbf{t}_3 \end{matrix}}\right] =
\left[{\begin{matrix}
\sigma _{11} & \sigma _{21} & \sigma _{31} \\
\sigma _{12} & \sigma _{22} & \sigma _{32} \\
\sigma _{13} & \sigma _{23} & \sigma _{33} 
\end{matrix}}\right] \cdot \left[{\begin{matrix} n_1 \\ n_2 \\ n_3 \end{matrix}}\right]
\quad \text{or} \quad
\mathbf{t} = \boldsymbol{\sigma} \cdot \mathbf{n}$$

This equation implies that the traction vector $\mathbf{t}$ at any point $P$ in a continuum associated with a plane given by unit normal vector $\mathbf{n}$, can be expressed as a function of the traction vectors on the planes perpendicular to the coordinate axes, i.e. stress tensor $\boldsymbol{\sigma}$.

```{code-cell} ipython3
S = stress.from_comp(xx=8, yy=6, zz=2)
S
```

```{code-cell} ipython3
n = fol(150, 60)  # normal of plane
t = S.cauchy(n)  # traction vector
print(f'Magnitude of normal stress on plane {n} is {abs(t.project(n))}')
print(f'Magnitude of shear stress on plane {n} is {abs(t.reject(n))}')
```

In APSG you can use the `stress_comp` method of the stress tensor:

```{code-cell} ipython3
sn, tau = S.stress_comp(n)
print(f'Magnitude of normal stress on plane {n} is {abs(sn)}')
print(f'Magnitude of shear stress on plane {n} is {abs(tau)}')
```

## Sign convention

```{figure} figures/konvence-tensors.png
:name: fig-sign-convention
:alt: Sign convention for stress components
:class: bg-primary mb-1
:width: 75%
:align: center

Sign convention for stress tensor components.
```

```{admonition} Sign convention for stress components
:class: tip
A stress component is positive if it acts in the positive direction of the coordinate axes, and if the plane where it acts has an outward normal vector pointing in the positive coordinate direction ({numref}`fig-sign-convention`).

Note that this rule by itself is convention-neutral — whether a positive component means *tension* or *compression* depends on which direction you take as "the" applied load. We follow the geosciences and rock-mechanics convention used by APSG: **compression is positive, tension is negative**. This is the opposite of the tension-positive convention common in engineering/continuum-mechanics texts, so take care when cross-referencing other sources.
```

## Mohr circle

The sign convention above tells us how to read a single traction vector. The **Mohr circle** is a graphical device that shows the normal and shear stress on *every* possible plane through a point in one picture, using the tensor's principal stresses. For a 2D stress state with principal stresses $\sigma_1 \ge \sigma_2$, the pair $(\sigma_n, \tau)$ for a plane whose normal makes angle $\theta$ with the $\sigma_1$ direction always lies on a circle:

$$\left(\sigma_n - \frac{\sigma_1+\sigma_2}{2}\right)^2 + \tau^2 = \left(\frac{\sigma_1-\sigma_2}{2}\right)^2$$

i.e. a circle centered at $\left(\frac{\sigma_1+\sigma_2}{2}, 0\right)$ with radius $\frac{\sigma_1-\sigma_2}{2}$. A key, easy-to-miss detail: rotating the physical plane by angle $\theta$ moves the point **twice as far** around the Mohr circle, i.e. by angle $2\theta$ — this "doubling" falls directly out of the $\cos(2\theta)$, $\sin(2\theta)$ terms that appear once you expand the tensor transformation rule (next chapter) in 2D.

We check both facts directly against a stress tensor, without asserting any formula blindly: sample many plane orientations, compute the true normal/shear stress from `.cauchy(n)` for each, and confirm every point lands on the predicted circle. (Note we build a fresh 2D tensor `S2d` here rather than reusing `S` from above, which by this point in the notebook has been reassigned to the 3D example in *Cauchy's stress theorem*.)

```{code-cell} ipython3
S2d = stress2([[8, -2],
               [-2,  5]])
s1, s2 = S2d.sigma1, S2d.sigma2
center, radius = (s1 + s2) / 2, (s1 - s2) / 2

thetas = np.linspace(0, 180, 91)
sn_vals, tau_vals = [], []
for theta in thetas:
    n = vec2(theta)
    t = S2d.cauchy(n)
    sn_vals.append(t.dot(n))     # signed normal stress
    tau_vals.append(t.cross(n))  # signed shear stress (2D perpendicular component)

fig, ax = plt.subplots()
circle = plt.Circle((center, 0), radius, fill=False, color='C1', label='predicted circle')
ax.add_patch(circle)
ax.plot(sn_vals, tau_vals, '.', color='C0', label='computed from S2d.cauchy(n)')
ax.axhline(0, color='gray', lw=0.5)
ax.set_xlabel(r'$\sigma_n$')
ax.set_ylabel(r'$\tau$')
ax.set_aspect(1)
ax.legend()
ax.set_title('Mohr circle for S2d');
```

Every computed point falls exactly on the predicted circle, and the point starting at $\theta=0$ reaches the opposite side of the circle at $\theta=90°$ — a $180°$ move on the circle for a $90°$ physical rotation, confirming the double-angle rule.

The *sign* of $\tau$ on a Mohr diagram is conventionally tied to a physical sense of shear (dextral vs. sinistral), depending on which way you view the plane — {numref}`fig-mohr-sense` shows the standard convention. We won't fix a specific pairing here (it depends on the viewing direction you choose), but the sign is not arbitrary: once a viewing convention is fixed, it consistently tracks one shear sense or the other.

```{figure} figures/konvence-mohr.png
:name: fig-mohr-sense
:alt: Mohr diagram shear sense convention
:class: bg-primary mb-1
:width: 90%
:align: center

Reading shear sense (dextral/sinistral) from the sign of $\tau$ on a Mohr diagram.
```

```{note}
The 2D construction above generalizes to 3D as **three** Mohr circles, one per pair of principal stresses ($\sigma_1,\sigma_2$), ($\sigma_2,\sigma_3$), ($\sigma_1,\sigma_3$); the admissible $(\sigma_n,\tau)$ region for *any* plane orientation is the area inside the outer circle and outside the two inner circles. We build this once we have a worked 3D example with known principal stresses in the next chapter.
```

