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

# Stress tensor

## Cauchy formula in 2D

```{image} figures/Cauchy_2D_triangle_v2.png
:alt: Cauchy in 2D
:class: bg-primary mb-1
:width: 50%
:align: center
```

The equilibrium of forces, i.e. Euler’s first law of motion (Newton’s second law of motion) for 2D triangle, gives:

$$\begin{aligned}
\boldsymbol{T}^{\left( n \right)} \cdot A &= \boldsymbol{T}^{\left( 1 \right)} \cdot A \cdot n_1 + \boldsymbol{T}^{\left( 2 \right)} \cdot A \cdot n_2\\
\boldsymbol{T}^{\left( n \right)} &= \boldsymbol{T}^{\left( 1 \right)} \cdot n_1
+ \boldsymbol{T}^{\left( 2 \right)} \cdot n_2
\end{aligned}$$

written by components we obtain

$$\begin{aligned}
T_1^{\left( n \right)} &= \sigma_{11} \cdot n_1 + \sigma_{21} \cdot n_2 \\
T_2^{\left( n \right)} &= \sigma_{12} \cdot n_1 + \sigma_{22} \cdot n_2 
\end{aligned}$$

or

$$\begin{bmatrix}
T_1^{\left( n \right)}\\ 
T_2^{\left( n \right)}
\end{bmatrix}
=
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
\boldsymbol{T}^{(\boldsymbol{n})} = \sigma _{ij} \cdot \bf{n}  \quad \text{where} \quad \sigma _{ij} =
\begin{bmatrix}
\boldsymbol{T}^{\left( 1 \right)}
&
\boldsymbol{T}^{\left( 2 \right)}
\end{bmatrix}$$

```{code-cell} ipython3
S = stress2([[-8,  2],
             [ 2, -5]])
n = vec2(60)  # unit length vector oriented 60 degrees from x
S.cauchy(n)  # traction vector
```

## Cauchy formula in 3D

The tetrahedron is formed by slicing the infinitesimal element along an arbitrary plane n. The stress vector on this plane is denoted by $T^{(n)}$. The stress vectors acting on the faces of the tetrahedron are denoted as $T^{(e_1)}$, $T^{(e_2)}$, and $T^{(e_3)}$, and are by definition the components $\sigma_{ij}$ of the stress tensor $\sigma$.

```{image} figures/Cauchy_tetrahedron.png
:alt: Cauchy in 3D
:class: bg-primary mb-1
:width: 60%
:align: center
```

The equilibrium of forces, i.e. Euler’s first law of motion (Newton’s second law of motion), gives:

$$\boldsymbol {T} ^{(\boldsymbol {n} )}\,dA-\boldsymbol {T} ^{(\boldsymbol {e} _{1})}\,dA_{1}-\boldsymbol {T} ^{(\boldsymbol {e} _{2})}\,dA_{2}-\boldsymbol {T} ^{(\boldsymbol {e} _{3})}\,dA_{3}=\rho \left({\frac {h}{3}}dA\right)\boldsymbol {a}$$

where the right-hand-side represents the product of the mass enclosed by the tetrahedron and its acceleration: $\rho$ is the density, $a$ is the acceleration, and $h$ is the height of the tetrahedron, considering the plane $n$ as the base.

The area of the faces of the tetrahedron perpendicular to the axes can be found by projecting dA into each face:

$$\begin{aligned}
dA_{1}&=\left(\boldsymbol {n} \cdot \boldsymbol {e} _{1}\right)dA=n_{1}\;dA\\
dA_{2}&=\left(\boldsymbol {n} \cdot \boldsymbol {e} _{2}\right)dA=n_{2}\;dA\\
dA_{3}&=\left(\boldsymbol {n} \cdot \boldsymbol {e} _{3}\right)dA=n_{3}\;dA
\end{aligned}$$

and then substituting into the equation to cancel out $dA$:

$$\boldsymbol {T} ^{(\boldsymbol {n} )}-\boldsymbol {T} ^{(\boldsymbol {e} _{1})}n_{1}-\boldsymbol {T} ^{(\boldsymbol {e} _{2})}n_{2}-\boldsymbol {T} ^{(\boldsymbol {e} _{3})}n_{3}=\rho \left({\frac {h}{3}}\right)\boldsymbol {a}$$

To consider the limiting case as the tetrahedron shrinks to a point, $h$ and the right-hand-side of the equation approaches 0, so:

$$\boldsymbol {T} ^{(\boldsymbol {n} )}=\boldsymbol {T} ^{(\boldsymbol {e} _{1})}n_{1}+\boldsymbol {T} ^{(\boldsymbol {e} _{2})}n_{2}+\boldsymbol {T} ^{(\boldsymbol {e} _{3})}n_{3}$$

## Cauchy stress tensor

**Tensors** are algebraic objects that describe linear relationship between vectors, scalars, or tensors. Here, any linear connection between two physical vector quantities is called a **tensor**, reflecting original use to describe the "tensions" in a material (Cauchy).

$$\left[{\begin{matrix} u_1 \\ u_2 \\ u_3 \end{matrix}}\right] = \left[{\begin{matrix}
    T_{11} & T_{21} & T_{31} \\
    T_{12} & T_{22} & T_{32} \\
    T_{13} & T_{23} & T_{33}
\end{matrix}}\right] \cdot \left[{\begin{matrix} v_1 \\ v_2 \\ v_3 \end{matrix}}\right]$$

$$\boldsymbol{u} = \boldsymbol{T} \cdot \boldsymbol{v}$$

In continuum mechanics, the Cauchy **stress tensor** $\boldsymbol\sigma$ is a second order tensor,
with nine components $\sigma_{ij}$, that completely define the state of stress at a point inside a material. The stress tensor is symmetric, so the number of independent stress components is equal to 6.

```{image} figures/Components_stress_tensor_cartesian.png
:alt: Stress tensor components
:class: bg-primary mb-1
:width: 60%
:align: center
```

$$\boldsymbol{\sigma} = \sigma_{ij} = \left[{\begin{matrix} \boldsymbol{T}^{(\boldsymbol{e}_1)} & \boldsymbol{T}^{(\boldsymbol{e}_2)} & \boldsymbol{T}^{(\boldsymbol{e}_3)} \end{matrix}}\right]  = \left[{\begin{matrix}
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

According to Cauchy’s fundamental theorem, also called Cauchy's stress theorem, merely by knowing the stress vectors on three mutually perpendicular planes,
the stress vector on any other plane passing through that point can be found through coordinate transformation equations.

$$\boldsymbol{T}^{(\boldsymbol{n})} \equiv \left[{\begin{matrix} {T_1}^{(\boldsymbol{n})} \\ {T_2}^{(\boldsymbol{n})} \\ {T_3}^{(\boldsymbol{n})} \end{matrix}}\right] = \left[{\begin{matrix}
\sigma _{11} & \sigma _{21} & \sigma _{31} \\
\sigma _{12} & \sigma _{22} & \sigma _{32} \\
\sigma _{13} & \sigma _{23} & \sigma _{33} 
\end{matrix}}\right] \cdot \left[{\begin{matrix} n_1 \\ n_2 \\ n_3 \end{matrix}}\right] \quad \text{or} \quad \boldsymbol{T}^{(\boldsymbol{n})} = \boldsymbol{\sigma} \cdot \bf{n}$$
This equation implies that the stress vector $\bf{T}^{\left( n \right)}$ at any point $P$ in a continuum associated with a plane with normal unit vector $\bf{n}$, can be expressed as a function of the stress vectors on the planes perpendicular to the coordinate axes, i.e. stress tensor $\boldsymbol{\sigma}$.

```{code-cell} ipython3
S = stress.from_comp(xx=-8, yy=-6, zz=-2)
S
```

```{code-cell} ipython3
n = fol(150, 60)  # normal of plane
T = S.cauchy(n)  # traction vector
print(f'Magnitude of normal stress on plane {n} is {abs(T.proj(n))}')
print(f'Magnitude of shear stress on plane {n} is {abs(T.reject(n))}')
```

## Sign convention

```{image} figures/konvence-tensors.png
:alt: Sign convention for stress components
:class: bg-primary mb-1
:width: 75%
:align: center
```

```{admonition} Sign convention for stress components
:class: tip
A stress component is positive if it acts in the positive direction of the coordinate axes, and if the plane where it acts has an outward normal vector pointing in the positive coordinate direction.
```

