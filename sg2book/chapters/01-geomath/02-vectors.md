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

# Vectors

Most structural geology students have learned the basics of vectors in their math courses.
We’ll first review those basic concepts and then put that knowledge to work because, as
mentioned in the last chapter, most linear features that we might wish to measure in
structural geology are vectors.

Because structural geometry is three dimensional, all of our vectors will have three
components. Each of the three numbers that define a vector refer to a specific coordinate axis:
for example, $u_2$ (or $u_y$) is the value of our vector, $\boldsymbol{u}$, projected onto the second axis of the
coordinate system, $X_2$ or $Y$. In a **NED** coordinate system, $u_2$ is the projection of $\boldsymbol{u}$ onto the
East axis, but in an **ENU** coordinate system, $u_2$ is u projected onto the North axis. Therefore,
the numbers that define a vector depend on the specific coordinate system. We write out
vector as:

$$\bf{u}=\begin{bmatrix} u_1 & u_2 & u_3 \end{bmatrix} = \begin{bmatrix} u_x & u_y & u_z \end{bmatrix}$$

```{code-cell} ipython3
u = vec(2, -1, 3)
v = vec(1, 2, -1)
```

## Vector magnitude

One of the most fundamental characteristics of a vector is its length or **magnitude**. Magnitude is a scalar quantity because it has no directional significance and furthermore, it is the same in all coordinate systems. The magnitude is given by:

$$\left | \boldsymbol{u} \right |=\sqrt{{u_1}^2+{u_2}^2+{u_3}^2}$$

In two dimensions, you can see that the magnitude is calculated from the Pythagorean theorem which gives the length of a hypotenuse as the square root of the sum of the squares of the two sides. The extension to three dimensions is straightforward.

```{image} figures/Vec_magnitude.png
:alt: Vector magnitude
:class: bg-primary mb-1
:width: 75%
:align: center
```

```{code-cell} ipython3
abs(u)
```

## Unit vector

```{image} figures/3D_Vector.png
:alt: 3D vector
:class: bg-primary mb-1
:width: 35%
:align: center
```

But, what if we don’t care about the magnitude? What if we are only interested in the orientation of our vector? It is convenient to represent direction by **unit vector**, i.e. the vector with length of one.
Any vector could be normalized to **unit vector** by dividing each of it's components by its magnitude.

$${\boldsymbol{\hat u}} = \frac{{\boldsymbol{u}}}{{\left\| {\boldsymbol{u}} \right\|}} = {\boldsymbol{i}}\frac{{{u_1}}}{{\left\| {\boldsymbol{u}} \right\|}} + {\boldsymbol{j}}\frac{{{u_2}}}{{\left\| {\boldsymbol{u}} \right\|}} + {\boldsymbol{k}}\frac{{{u_3}}}{{\left\| {\boldsymbol{u}} \right\|}}$$

The projection of a unit vector onto a coordinate axis is just equal to the cosine of the angle that the vector makes with that axis.

Thus, the components of a unit vector are:

$${\boldsymbol{\hat u}} = \begin{bmatrix} \cos (\alpha ) & \cos (\beta ) & \cos (\gamma ) \end{bmatrix}$$

```{code-cell} ipython3
u.normalized()
```

## Vector addition, subtraction and scalar multiplication

Vector addition or subtraction from another is no more complicated than adding, or subtracting, the individual components:

$$\bf{u+v}=\begin{bmatrix} u_1 & u_2 & u_3 \end{bmatrix}+\begin{bmatrix} v_1 & v_2 & v_3\end{bmatrix}=\begin{bmatrix} (u_1+v_1) & (u_2+v_2) & (u_3+v_3) \end{bmatrix}$$

And likewise

$$\bf{u-v}=\begin{bmatrix} u_1 & u_2 & u_3 \end{bmatrix}-\begin{bmatrix} v_1 & v_2 & v_3\end{bmatrix}=\begin{bmatrix} (u_1-v_1) & (u_2-v_2) & (u_3-v_3) \end{bmatrix}$$

```{image} figures/Vec_add_sub.png
:alt: Vector addition and subtraction
:class: bg-primary mb-1
:width: 75%
:align: center
```

```{code-cell} ipython3
u + v
```

## Dot product

**Dot product** or **scalar product** is an algebraic operation between two vectors. The dot product may be defined algebraically or geometrically. Conceptually, **dot product** can be thought of as multiplying the length of one vector by the *component of the other vector* which is parallel to the first:

Geometrically, it is the product of the magnitudes of the two vectors and the cosine of the angle between them.

$$\boldsymbol u\cdot\boldsymbol v = \|\boldsymbol u\|\,\|\boldsymbol v\|\cos\theta$$

Algebraically, it is the sum of the products of the corresponding entries of the two sequences of numbers.

$$\boldsymbol{u}\cdot \boldsymbol{v} = \sum_{i=1}^n u_iv_i = u_1v_1 + u_2v_2 + \cdots + u_nv_n$$

Two vectors are orthogonal if the dot product of those two vectors is equal to zero.

```{code-cell} ipython3
u.dot(v)
```

## Vector projection

The **vector projection** of a vector $\boldsymbol{u}$ on a nonzero vector $\boldsymbol{v}$ (also known as the *vector component* of $\boldsymbol{u}$ in the direction of $\boldsymbol{u}$) is the orthogonal projection of $\boldsymbol{u}$ onto a straight line parallel to $\boldsymbol{v}$. It is a vector parallel to $\boldsymbol{v}$, defined as:

```{image} figures/Dot_Product.png
:alt: Vector projection
:class: bg-primary mb-1
:width: 35%
:align: center
```

$$\boldsymbol{u}_1 = u_1\boldsymbol{\hat v} = u_1\frac {\boldsymbol{v}} {\|\boldsymbol{v}\|},$$

where $u_1$ is **scalar projection**. Using definition of cosine in right-angled triangle and **dot product** definition, we can write:

$$u_1 = \|\boldsymbol{u}\| \cos \theta = \frac {\boldsymbol{u} \cdot \boldsymbol{v}} {\|\boldsymbol{v}\| }$$

Consequently,

$$\boldsymbol{u}_1 = u_1\frac {\boldsymbol{v}} {\|\boldsymbol{v}\|} = \frac {\boldsymbol{u} \cdot \boldsymbol{v}} {\|\boldsymbol{v}\| } \frac {\boldsymbol{v}} {\|\boldsymbol{v}\|} = \frac {\boldsymbol{u} \cdot \boldsymbol{v}} {\|\boldsymbol{v}\|^2}{\boldsymbol{v}} = \frac {\boldsymbol{u} \cdot \boldsymbol{v}} {\boldsymbol{v} \cdot \boldsymbol{v}}{\boldsymbol{v}}$$

```{code-cell} ipython3
u.proj(v)
```

## Cross product

The **cross product** $\boldsymbol{u} \times \boldsymbol{v}$ is defined as a vector that is perpendicular to both $\boldsymbol{u}$ and $\boldsymbol{v}$, with a direction given by the right-hand rule and a magnitude equal to the area of the parallelogram that the vectors span.

```{image} figures/Cross_product_vector.png
:alt: Cross product
:class: bg-primary mb-1
:width: 25%
:align: center
```

The cross product is defined by the formula

$$\boldsymbol{u} \times \boldsymbol{v} = \left\| \boldsymbol{u} \right\| \left\| \boldsymbol{v} \right\| \sin \theta \ \boldsymbol{n}$$

or in **matrix notation**

$$\boldsymbol{u\times v}=\begin{vmatrix} \boldsymbol{i}&\boldsymbol{j}&\boldsymbol{k}\\ u_1&u_2&u_3\\ v_1&v_2&v_3\\ \end{vmatrix}$$

The cross product is *anticommutative*

$$\boldsymbol{u} \times \boldsymbol{v} = -(\boldsymbol{v} \times \boldsymbol{u})$$

```{code-cell} ipython3
u.cross(v)
```

## Vectors in geology

The compass measurements are commonly in **spherical coordinates** i.e. *trend*
or *dip direction* ($\theta$) and *plunge* or *dip* ($\varphi$). *Plunge* is used
to describe the tilt of lines, the word *dip* being reserved for planes.

```{image} figures/Geo_coordinates_en.png
:alt: Coordinate transformation in geology
:class: bg-primary mb-1
:width: 50%
:align: center
```

The linear feature has (a) the plunge direction or trend and (b) the angle of plunge.
The plunge direction is the direction towards which the line is tilted.
The angle of plunge is the amount of tilt; it is the angle, measured in the
vertical plane, that the plunging line makes with the horizontal. The angle of
plunge of a horizontal line is 0° and the angle of plunge of a vertical line is 90°. 

Linear features (lines) are represented as vector:

$$\begin{aligned}
u_1&=&\cos(\theta)\cos(\varphi)\\
u_2&=&\sin(\theta)\cos(\varphi)\\
u_3&=&\sin(\varphi)
\end{aligned}$$

Planar fetaure that is not horizontal is said to dip. There are two aspects to
the dip of a plane (a) the direction of dip, which is the compass direction
towards which the plane slopes; and (b) the angle of dip, which is the angle
that the plane makes with a horizontal plane. The direction of dip can be
visualized as the direction in which water would flow if poured onto the
plane. The angle of dip is an angle between 0° (for horizontal planes) and 90°
(for vertical planes).

Planar features (planes) are represented by normal vector:

$$\begin{aligned}
p_1&=&-\cos(\theta)\sin(\varphi)\\
p_2&=&-\sin(\theta)\sin(\varphi)\\
p_3&=&\cos(\varphi)
\end{aligned}$$

With APSG the 3D vector could be defined using trend and plunge notation:

```{code-cell} ipython3
v = vec(130, 20)
```

while linear and planar features could be defined as lineation or foliation

```{code-cell} ipython3
l = lin(120, 40)
f = fol(210, 30)
quicknet(f, l, fol_as_pole=False)
```

