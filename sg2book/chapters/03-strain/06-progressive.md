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
from scipy.linalg import expm
from apsg import *
from strain2d import plot_defgrad
```

# Progressive deformation

Everything so far has compared a single **reference** configuration to a single **deformed** configuration through one deformation gradient $\boldsymbol{F}$. Real geological deformation, however, accumulates continuously over time. This chapter follows {cite}`means1976` and {cite}`ramsayhuber1983` in describing that continuous accumulation, and revisits **coaxial** vs. **non-coaxial** deformation (introduced in *Superposition of deformation*) and **pure shear** vs. **simple shear** (introduced in *Deformation gradient decompositions*) as end-members of a single continuous quantity: the **kinematic vorticity number**.

## Velocity gradient

If $\mathbf{x}(t)$ is the position of a material point at time $t$, its velocity is $\dot{\mathbf{x}}=d\mathbf{x}/dt$. The **velocity gradient** $\boldsymbol{L}$ relates velocity to the current (not reference) configuration:

$$\boldsymbol{L} = \dot{\boldsymbol{F}}\boldsymbol{F}^{-1}$$

Unlike $\boldsymbol{F}$, which compares two finite states, $\boldsymbol{L}$ describes the *instantaneous rate* of deformation at a single point in time. Like any second-order tensor, it splits uniquely into symmetric and antisymmetric parts:

$$\boldsymbol{L} = \underbrace{\tfrac{1}{2}\left(\boldsymbol{L}+\boldsymbol{L}^T\right)}_{\boldsymbol{D}\ \text{(strain rate)}} + \underbrace{\tfrac{1}{2}\left(\boldsymbol{L}-\boldsymbol{L}^T\right)}_{\boldsymbol{W}\ \text{(spin / vorticity)}}$$

$\boldsymbol{D}$ describes the instantaneous rate of stretching (its eigenvectors are the instantaneous stretching axes); $\boldsymbol{W}$ describes the instantaneous rigid rotation rate of those axes — exactly the same symmetric/antisymmetric split as $\boldsymbol{\nabla u} = \boldsymbol{\varepsilon}+\boldsymbol{\omega}$ from the infinitesimal-strain discussion two chapters back, but now applied to a *rate* rather than to a single finite displacement.

## Steady flow and the deformation path

If $\boldsymbol{L}$ is constant in time (a **steady** progressive deformation), the finite deformation gradient accumulated after time $t$ is given by the matrix exponential:

$$\boldsymbol{F}(t) = \exp(\boldsymbol{L}t)$$

This lets us recover the pure-shear and simple-shear examples from the previous chapter as the two simplest possible steady flows:

```{code-cell} ipython3
# steady pure shear: L diagonal (no spin) — constant stretching rate
k = np.log(2.0)
L_pure = np.array([[k, 0], [0, -k]])
F_pure = expm(L_pure * 1.0)
print(F_pure)  # matches U = [[2, 0], [0, 0.5]] from the polar-decomposition chapter
```

```{code-cell} ipython3
# steady simple shear: L off-diagonal only — constant shear rate
L_simple = np.array([[0, 1], [0, 0]])
F_simple = expm(L_simple * 1.0)
print(F_simple)  # matches Fs = [[1, 1], [0, 1]] from the polar-decomposition chapter
```

Integrating a *constant* stretching-only $\boldsymbol{L}$ for one time unit reproduces exactly the stretch tensor $\boldsymbol{U}$ used earlier, and integrating a constant shear-only $\boldsymbol{L}$ reproduces exactly the simple-shear $\boldsymbol{F}_s$ used earlier — pure shear and simple shear are not just two isolated matrix examples, they are the finite results of the two simplest possible steady flows.

## Kinematic vorticity number

The ratio of spin to stretching rate along the principal stretching axes is the **kinematic vorticity number** {cite}`means1976`:

$$W_k = \frac{|\omega|}{\dot{e}}$$

where $\omega$ is the (single, in 2D) independent component of $\boldsymbol{W}$ and $\dot{e}$ is the magnitude of the principal strain rate (the largest eigenvalue of $\boldsymbol{D}$, for an isochoric flow with eigenvalues $\pm\dot{e}$).

```{code-cell} ipython3
def kinematic_vorticity(L):
    D = 0.5 * (L + L.T)
    W = 0.5 * (L - L.T)
    edot = np.linalg.eigvalsh(D).max()
    omega = abs(W[1, 0])
    return omega / edot if edot != 0 else np.inf

print('Wk (pure shear)  =', kinematic_vorticity(L_pure))
print('Wk (simple shear) =', kinematic_vorticity(L_simple))
```

```{admonition} Wk as a single dial between pure and simple shear
:class: tip
- $W_k=0$: no spin at all — the stretching axes never rotate relative to material lines. This is exactly the definition of a **coaxial**, **pure-shear**-type flow.
- $W_k=1$: **simple shear** — material lines parallel to the shear plane are not rotated, but the principal stretching axes rotate continuously relative to the material. This is the standard reference case for **non-coaxial** flow.
- $0<W_k<1$: **general non-coaxial flow** (sometimes called *sub-simple shear*), a combination of pure and simple shear components with the same principal directions.
```

## Visualizing a progressive, non-coaxial deformation path

A single snapshot $\boldsymbol{F}$ cannot show *how* the ellipse got there. For a steady flow, we can evaluate $\boldsymbol{F}(t)=\exp(\boldsymbol{L}t)$ at a sequence of times and watch the strain ellipse evolve — this is exactly what a progressive, non-coaxial deformation looks like in the field over successive increments.

```{code-cell} ipython3
# general non-coaxial flow with Wk = 0.5 (between pure and simple shear)
e, omega = 0.5, 0.25
L = np.array([[e, omega], [-omega, -e]])
print('Wk =', kinematic_vorticity(L))

fig, axs = plt.subplots(1, 4, figsize=(12, 3))
for ax, t in zip(axs, [0.25, 0.5, 0.75, 1.0]):
    F_t = defgrad2(expm(L * t).tolist())
    plot_defgrad(F_t, ax=ax, title=f't={t:g}')
```

Unlike the pure-shear case, the ellipse's long axis does not stay fixed relative to the reference frame as it grows — it keeps rotating throughout the deformation history, the visual signature of $0<W_k<1$.

```{note}
Composing many small increments of a constant $\boldsymbol{L}$ — $\boldsymbol{F}(t+\Delta t)\approx(\boldsymbol{I}+\boldsymbol{L}\Delta t)\cdot\boldsymbol{F}(t)$, applied repeatedly — converges to $\exp(\boldsymbol{L}t)$ as $\Delta t\to 0$. This is the same non-commuting composition $\boldsymbol{F}_2\cdot\boldsymbol{F}_1\ne\boldsymbol{F}_1\cdot\boldsymbol{F}_2$ from *Superposition of deformation*, taken to its continuous limit: a non-coaxial progressive deformation is, fundamentally, an infinite sequence of tiny non-commuting increments.
```
