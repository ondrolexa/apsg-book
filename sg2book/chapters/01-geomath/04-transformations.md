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

# Transformations

Informally, *function* is a rule that accepts inputs and produces outputs. For instance, $f(x)=x^2$ is a function that accepts one number $x$ as its input, and outputs the square of that number: $f(2)=4$.

Consider the matrix equation:

$$\mathbf{v} = \mathbf{A}\mathbf{u}$$

If we vary $\mathbf{u}$, then $\mathbf{v}$ will also vary; in this way, we think of $\mathbf{A}$ as *function* with independent variable $\mathbf{u}$ (which is a vector in $\mathbb{R}^n$) and dependent variable $\mathbf{v}$ (which is a vector in $\mathbb{R}^m$).

At this point it is convenient to fix our ideas and terminology regarding *functions*, which we will call **transformations**. This allows us to systematize our discussion of **matrices as functions**.

```{code-cell} ipython3
u = vec(2, -1, 3)
F = matrix([[2, 1, 0], [0, 1, 0], [0, 0, 0.5]])
u.transform(F)
```

## Matrices as functions - transformation matrix

If $\mathbf{A}$ has $n$ columns, then it only makes sense to multiply $\mathbf{A}$ by vectors with $n$ entries. This is why the domain of transformation is $\mathbb{R}^n$. If $\mathbf{A}$ has $m$ rows, then $\mathbf{A}\mathbf{u}$ has $m$ entries for any vector $\mathbf{u}$ in $\mathbb{R}^n$; this is why the codomain of transformation is $\mathbb{R}^m$.

Suppose that $\mathbf{A}$ has columns $v_1, v_1, \dots v_n$,. If we multiply $\mathbf{A}$ by a general vector $\mathbf{u}$, we get:

$$\mathbf{A}\mathbf{u}=\begin{pmatrix}| & | &  & | \\ v_1 & v_2 & \dots & v_3\\ | & | &  & | \end{pmatrix}\begin{pmatrix}u_1\\ u_2\\ \vdots \\ u_n\end{pmatrix}= u_1v_1 + u_2v_2 + \dots + u_nv_n$$

This is just a general linear combination of $v_1, v_1, \dots v_n$. Therefore, the outputs of matrix transformations are exactly the linear combinations of the columns of *transformation matrix* $\mathbf{A}$.

```{admonition} Useful identities
:class: tip
Following identities are often used in linear algebra:

$$\begin{aligned}
(\mathbf{A}\mathbf{B})^T &= \mathbf{B}^T\mathbf{A}^T\\
(\mathbf{A}\mathbf{B})^{-1} &= \mathbf{B}^{-1}\mathbf{A}^{-1}\\
(\mathbf{A}^T)^{-1} &= (\mathbf{A}^{-1})^T = \mathbf{A}^{-T}\\
((\mathbf{A}\mathbf{B})^{-1})^T &= ((\mathbf{A}\mathbf{B})^T)^{-1} = (\mathbf{B}^{-1}\mathbf{A}^{-1})^T = \mathbf{A}^{-T}\mathbf{B}^{-T}
\end{aligned}$$
```