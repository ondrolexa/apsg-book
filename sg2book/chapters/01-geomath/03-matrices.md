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

# Matrices

A matrix, like a vector, is also a collection of numbers. The difference is that a matrix is a table of numbers rather than a list. Many of the same rules we just outlined for vectors above apply equally well to matrices. In fact, you can think of vectors as matrices that happen to only have one column or one row.

The **dimensions of a matrix** tells its size: the number of rows and columns of the matrix, in that order.

Since matrix $\boldsymbol{A}$ has two rows and three columns, we write its dimensions as $2\times 3$, pronounced "two by three". In contrast, matrix $\boldsymbol{B}$ has three rows and two columns, so it is a $3\times 2$ matrix.

$$\boldsymbol{A} = \begin{bmatrix}
  -2 & 5 & 6\\ 
  5 & 2 & 7
\end{bmatrix}$$

$$\boldsymbol{B} = \begin{bmatrix}
  -8 & -4\\ 
  23 & 12\\ 
  18 & 10
\end{bmatrix}$$

```{code-cell} ipython3
A = matrix2([[-2, 5], [5, 2]])
B = matrix2([[4, -1], [-1, 3]])
```

## Matrix addition, substraction and multiplication

First, let’s consider matrix addition and subtraction. This part is uncomplicated. You can add and subtract matrices the same way you add vectors – element by element:

$$\boldsymbol{A}+\boldsymbol{B}=\boldsymbol{C}$$

$$\begin{bmatrix} a_{11} & a_{12}\\a_{21} & a_{22}\end{bmatrix} + \begin{bmatrix} b_{11} & b_{12}\\b_{21} & b_{22}\end{bmatrix} = \begin{bmatrix} a_{11}+b_{11} & a_{12}+b_{12}\\a_{21}+b_{21} & a_{22}+b_{22}\end{bmatrix}$$

```{code-cell} ipython3
A + B
```

Matrix multiplication gets a bit more complicated, since multiple elements in the first matrix interact with multiple elements in the second to generate each element in the product matrix. This means that matrix multiplication can be a tedious task to carry out by hand, and can be time consuming on a computer for very large matrices.

$$\boldsymbol{A}\boldsymbol{B}=\boldsymbol{C}$$

$$\begin{bmatrix} a_{11} & a_{12}\\a_{21} & a_{22}\end{bmatrix} \begin{bmatrix} b_{11} & b_{12}\\b_{21} & b_{22}\end{bmatrix} = \begin{bmatrix} a_{11}b_{11}+a_{12}b_{21} & a_{11}b_{12}+a_{12}b_{22}\\a_{21}b_{21}+a_{12}b_{21} & a_{21}b_{22}+a_{12}b_{22}\end{bmatrix}$$

```{code-cell} ipython3
A @ B
```

