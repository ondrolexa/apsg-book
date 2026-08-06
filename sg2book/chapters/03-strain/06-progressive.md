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
from scipy.linalg import expm, logm
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

## From 2D to 3D: APSG's velocity and deformation gradients

Everything above used plain $2\times2$ arrays and `scipy.linalg.expm` directly, so that the algebra stayed visible. APSG's `velgrad`/`defgrad` classes wrap exactly the same $\boldsymbol{F}=\exp(\boldsymbol{L}t)$ relationship in three dimensions, and are the more practical tool once $\boldsymbol{L}$ is no longer a hand-picked $2\times2$ matrix — they already appear elsewhere in this book, e.g. `ellipse.from_defgrad` in *Deformation gradient decompositions*. {cite}`davistitus2011` give a comprehensive, coordinate-free review of this exponential/logarithm framework for homogeneous steady deformation, which the rest of this chapter follows.

```{code-cell} ipython3
# the same steady pure shear as above, now as a full 3x3 tensor (zz chosen to preserve volume)
L3 = velgrad.from_comp(xx=k, yy=-k, zz=0)
F3 = L3.defgrad(time=1)
print(F3)
print(F3.velgrad())  # round-trip: recovers L3 exactly
```

## Non-uniqueness of the matrix logarithm

The round trip above ran in a single, reassuring direction: build $\boldsymbol{L}$, exponentiate, get $\boldsymbol{F}$. The reverse direction — recovering $\boldsymbol{L}$ from a known $\boldsymbol{F}$ via $\boldsymbol{L}=\ln\boldsymbol{F}$ — is the one geologists actually need, since field measurements record finite states, not velocity fields, and it is not always safe. {cite}`davistitus2011` show that two different velocity gradients can integrate to the very same finite deformation: one purely coaxial, the other identical but with an extra rigid rotation added that never shows up in the end state, because a $2\pi$ rotation returns every material line to where it started.

```{admonition} Only the principal logarithm is recovered
:class: tip
`scipy.linalg.logm` (and APSG's `.velgrad()`) always return the **principal** logarithm of $\boldsymbol{F}$ — the one with the least rotation. If the true deformation path underwent one or more extra full rotations along the way, $\ln\boldsymbol{F}$ recovers "the simplest explanation" of the finite state, not necessarily the true kinematic history. Rotation, in general, cannot be recovered from the finite state alone.
```

```{code-cell} ipython3
# Two velocity gradients producing the SAME finite deformation (Davis and Titus, 2011, Fig. 4)
lam, n_turns = 0.5, 1

L0 = velgrad.from_comp(xx=-lam, yy=-lam, zz=2 * lam)  # purely coaxial
L1 = velgrad([[-lam, -2 * np.pi * n_turns, 0],
              [2 * np.pi * n_turns, -lam, 0],
              [0, 0, 2 * lam]])                        # L0 plus a rotation

F0 = L0.defgrad(time=1)
F1 = L1.defgrad(time=1)
print(np.allclose(F0, F1))  # same F, from two different L
```

```{code-cell} ipython3
print(F1.velgrad())  # recovers L0, the principal logarithm -- not L1
```

## Simultaneous superposition: the ln-sum-exp method

*Superposition of deformation* showed that composing two finite deformations sequentially, $\boldsymbol{F}=\boldsymbol{F}_2\boldsymbol{F}_1$, is order-dependent: applying $\boldsymbol{F}_1$ then $\boldsymbol{F}_2$ is generally not the same as applying $\boldsymbol{F}_2$ then $\boldsymbol{F}_1$. But some field situations call for two deformation processes acting *at the same time* rather than one after another — a shear zone widening through boundary-parallel slip while, simultaneously, the whole zone shortens to conserve volume, for example. If two processes proceed simultaneously, their velocity fields simply add: $\boldsymbol{L}=\boldsymbol{L}_1+\boldsymbol{L}_2$. Because $\boldsymbol{L}_i=\ln\boldsymbol{F}_i$ for each component, the finite deformation that results from running $\boldsymbol{L}_1$ and $\boldsymbol{L}_2$ together for one time unit is

$$\boldsymbol{F} = \exp(\ln\boldsymbol{F}_1+\ln\boldsymbol{F}_2)$$

— {cite}`davistitus2011`'s "ln-sum-exp" method, justified by the Trotter product formula linking this closed form to the limit of ever-finer alternating increments of $\boldsymbol{F}_1$ and $\boldsymbol{F}_2$. It is a genuinely different operation from the sequential $\boldsymbol{F}_2\boldsymbol{F}_1$ of the previous chapter, not merely a reformulation of it.

```{code-cell} ipython3
# F1: boundary-preserving transpression (shear + along-strike stretch, no volume correction)
# F2: a vertical stretch that exactly restores the volume F1 loses
v1, v2, v3 = 0.4, 0.3, 0.1

F1 = defgrad([[1, v1, 0], [0, 1 + v2, 0], [0, v3, 1]])
F2 = defgrad([[1, 0, 0], [0, 1, 0], [0, 0, 1 / (1 + v2)]])
print('det F1 =', F1.det, ' det F2 =', F2.det)
```

```{code-cell} ipython3
# run F1 and F2 simultaneously (not sequentially): sum their logs, then exponentiate
L = velgrad(np.asarray(F1.velgrad()) + np.asarray(F2.velgrad()))
F = L.defgrad()
print(F)
print('det F =', F.det)  # volume is exactly restored
```

```{code-cell} ipython3
print(np.allclose(F, F2 @ F1))  # NOT the same as sequential composition
```

## Application: inclined transpression and fabric development

Many natural shear zones show *oblique* lineations — neither strike-parallel nor dip-parallel — that classical monoclinic transpression models (a single simple shear plus a coaxial component sharing its principal axes) cannot produce. {cite}`jonesholdsworth1998` proposed an *inclined* transpression model with two independent simple-shear components on perpendicular planes plus a coaxial stretch, giving the deformation a genuinely triclinic symmetry. {cite}`davistitus2011` rebuild this model with the ln-sum-exp method above: summing the principal logarithms of two orthogonal simple shears (rates $\gamma_{xy}$, $\gamma_{zy}$) and a volume-preserving coaxial stretch ($\alpha_z$) gives, in one step, the velocity gradient

$$\boldsymbol{L}=\begin{bmatrix}0&\gamma_{xy}&0\\0&-\ln\alpha_z&0\\0&\gamma_{zy}&\ln\alpha_z\end{bmatrix}$$

```{code-cell} ipython3
gxy, gzy, az = 0.3, 0.2, 1.5  # two shear rates and a coaxial stretch factor

L = velgrad.from_comp(xy=gxy, yy=-np.log(az), zy=gzy, zz=np.log(az))
F = L.defgrad(time=1)
print(F)
print('det F =', F.det)  # volume-preserving, by construction
```

The kinematic vorticity number defined earlier in this chapter, $W_k=|\omega|/\dot e$, was written for a 2D flow with a single independent spin component. In 3D, {cite}`davistitus2011` give the coordinate-free definition

$$W_k=\frac{|\boldsymbol{W}|}{|\boldsymbol{S}|}$$

where $|\cdot|$ is the Frobenius norm and $\boldsymbol{S}=\tfrac12(\boldsymbol{L}+\boldsymbol{L}^T)$, $\boldsymbol{W}=\tfrac12(\boldsymbol{L}-\boldsymbol{L}^T)$ are exactly the $D$/$W$ split from the start of this chapter — the same quantity, generalized beyond two dimensions.

```{code-cell} ipython3
S, W = L.rate(), L.spin()
Wk = np.linalg.norm(np.asarray(W)) / np.linalg.norm(np.asarray(S))
print('Wk =', Wk)  # between pure shear (0) and simple shear (1): general non-coaxial transpression
```

A transpressional flow like this one leaves a distinctive fingerprint on a rock's fabric. Starting from an initially near-random distribution of orientations (e.g. detrital grains with no preferred orientation) and tracking the orientation tensor as the flow accumulates traces out a path on a fabric-shape diagram; here we use the Hsu diagram (octahedral shear strain against Lode's parameter), which stays legible even though this flow, unlike the field examples in earlier chapters, is not plane strain.

```{code-cell} ipython3
np.random.seed(1)
v = vecset.random(300, name='initial fabric')  # uniform, unstrained starting distribution

# orientation tensor evolution along the deformation path F(t), t in [0, 1]
path = ortensorset([ortensor.from_features(v.transform(L.defgrad(t)))
                     for t in np.linspace(0, 1, 11)])

h = HsuPlot()
h.path(path, marker='.')
h.point(path[-1], color='r', label='final fabric')
h.show()
```

As deformation accumulates, the initially weak, near-random fabric strengthens steadily and its shape drifts toward the oblate/flattening field (Woodcock's $K$ falls from just above 1 toward well below it). The shear components alone would drive the fabric toward a girdle, and the coaxial vertical stretch alone would drive it toward a point maximum; in this particular transpression, the girdle tendency wins — precisely the kind of coaxial/non-coaxial interplay $W_k$ was built to quantify.
