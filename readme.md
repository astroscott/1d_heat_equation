## 1D Heat Transfer In Aluminum Beam

The state, ie: temperature at some position $\mathbf{x}$ at some time $t$.

$$ u(\mathbf{x},t) $$ 

The state derivative is the heat equation:

$$\frac{\partial u}{\partial t} = \alpha \nabla^2 u $$

where

$u$ = temperature

$\alpha$ = thermal diffusivity

$t$ = time

Considering only 1-dimension:

$$ \nabla^2 u = \frac{\partial^2 u}{\partial x^2} = \frac{u_{i+1}^m - 2u_i^m + u_{i-1}^m}{\Delta x^2}$$

where $m$ represents a given time, and $i$ represents spatial discretization. All together, the discrete heat equation is:

$$ \therefore \frac{\partial u}{\partial t} = \alpha \left( \frac{u_{i+1}^m - 2u_i^m + u_{i-1}^m}{\Delta x^2} \right) $$

The inital condition and boundary conditions are defined as follows:

IC: $ u(x,0) = 273.15 K $ for $ 0 \lt x \lt L $, where L is the total length of the rod

BC1: $ u(0, t) = 500 K $ for all t

BC2: $ u(L, t) = 500 K $ for all t

The state derivative can be expressed in matrix form:

```math
\therefore \frac{\partial \mathbf{u}}{\partial t} = \frac{\alpha}{\Delta x^2}

\begin{bmatrix}
-2 & 1 & 0 & \dots & 0\\
1 & -2 & 1 & \ddots & 0\\
0 & 1 & -2 & \ddots & 0\\
\vdots & \ddots & \ddots & \ddots & 1 \\
0 & 0 & 0 & 1 & -2 

\end{bmatrix}
\begin{bmatrix}
u_0 \\
u_1 \\
u_2 \\
\vdots \\
u_n
\end{bmatrix}
```

and implemented using scipy.integrate.solve_ivp(), then plotted in matplotlib.