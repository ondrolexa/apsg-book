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

# Coordinate Systems

When we use a compass to measure features of interest to the structural geologist, we are
implicitly using a spherical coordinate system defined by the rotation axis and surface of the
Earth. Even if the Earth is almost spherical body, it is commonly more convenient to do
structural geology calculations in rectangular Cartesian coordinates.

The axes of any graph have positive and negative directions and most of us are used to seeing
two dimensional graphs where the horizontal axis X is positive to the right and the vertical Y
axis is positive upwards. But, what happens if we add a third axis? How do we determine the
positive direction for that axis?

The three axes of a Cartesian coordinate system are commonly referred to as $X$, $Y$, and $Z$. We
start out the same way here because that is most familiar, but when we start talking about
coordinate systems more formally, we will switch to using $X_1$, $X_2$, and $X_3$, which are more
convenient for numerical calculations in a computer.

## Coordinate systems in geology

The angles are measured positive in a direction from the $X$ axis towards the $Y$ axis, i.e., just
the opposite of how angles are measured on a compass rose (clockwise from the top).

Convention suggests that we should follow a right-handed naming convention:

Right-hand rule indicates the direction of the coordinate axes. When you hold the thumb, index finger, and middle finger of your right hand so that they form three right angles, then the thumb symbolizes the $x$ axis, the index finger the $y$ axis, and the middle finger the $z$ axis.

```{image} figures/Right_hand_rule.png
:alt: Vector projection
:class: bg-primary mb-1
:width: 50%
:align: center
```

If you hold you hand so that your thumb points in the positive direction of the first axis, your fingers should
curl from the positive direction of the second axis toward the positive direction of the third.

```{image} figures/Left-Right.png
:alt: Left-right-hand
:class: bg-primary mb-1
:width: 75%
:align: center
```

The two systems are commonly adopted in geology, **East-North-Up (ENU)** and **North-East-
Down (NED)** coordinate system.

```{image} figures/Geocoord_sytems.png
:alt: Geo-coord-systems
:class: bg-primary mb-1
:width: 75%
:align: center
```

Of course it is easy to change between the two coordinate systems: A point that has coordinates
(100, 200, 350) in an ENU system would have coordinates (200, 100, and -350) in a NED system.
