@[TOC](Contents)

<script type="text/javascript"
src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

## Vector Calculus and Partial Differential Equations

### 1.1 Big Picture Overview

Conservation Laws
* Mass
* Momentum
* Energy

Vector Calc
* $\nabla $   grad
* $\nabla \cdot$  div
* $\nabla \times$ curl

### 1.2 Div, Grad, and Curl: Vector Calculus Building Blocks for PDEs [Divergence, Gradient, and Curl]

$$\vec\nabla = \frac{\partial}{\partial x} \vec i + \frac{\partial}{\partial y} \vec j + \frac{\partial}{\partial z}\vec k \tag{1}$$

or

$$\vec \nabla = \left[\begin{array}{} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y} \\ \frac{\partial}{\partial z}  \end{array} \right]$$

is the gradient operator (grad operation),
- "Grad" gradient operator takes a scalar $f \to$ a vector field;
- "Div" divergence operatoer takes a vector $\vec f \to$ scalar field;
$$\nabla \cdot \vec f = \left[\begin{array}{} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y} \\ \frac{\partial}{\partial z}  \end{array} \right] \cdot \left[\begin{array}{} \vec i f_1 \\ \vec j f_2 \\ \vec k f_3 \end{array}\right] =  \frac{\partial f_1}{\partial x} + \frac{\partial f_2}{\partial y} + \frac{\partial f_3}{\partial z} $$
- "Curl" curl operator takes a vector $\vec f \to $ a vector field.
$$\nabla \times \vec f = \left[\begin{array}{} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y} \\ \frac{\partial}{\partial z}  \end{array} \right] \times \left[\begin{array}{} \vec i f_1 \\ \vec j f_2 \\ \vec k f_3 \end{array}\right] =  \left[ \begin{array}{ccc} \vec i & \vec j & \vec k \\ \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z} \\ f_1 & f_2 & f_3  \end{array} \right] = \vec i \left(\frac{\partial f_3}{\partial y}-\frac{\partial f_2}{\partial z}\right) - \vec j \left(\frac{\partial f_3}{\partial x}-\frac{\partial f_1}{\partial z}\right) \vec k \left(\frac{\partial f_2}{\partial x}-\frac{\partial f_1}{\partial y}\right) \tag{2}$$


## 2. The Gradient Operator in Vector Calculus

Directions of Fastest Change & the Directional Derivative

## 3. The Divergence of a Vector Field

Sources and Sinks

$$\nabla \cdot \vec f = \left[\begin{array}{} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y}  \end{array} \right] \cdot \left[\begin{array}{} \vec i f_1 \\ \vec j f_2 \end{array}\right] =  \frac{\partial f_1}{\partial x} + \frac{\partial f_2}{\partial y} $$

### 3.1 Divergence of gradient is the Laplacian

$$\vec \nabla\cdot \nabla f = \left[\begin{array}{} \frac{\partial}{\partial x} \\ \frac{\partial}{\partial y}  \end{array} \right] \cdot \left[\begin{array}{} \frac{\partial f}{\partial y} \\ \frac{\partial f}{\partial y} \end{array}\right] = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} = \nabla^2 f$$


## 4. The Curl of a Vector Field

Measuring Rotation

### 4.1 Curl(Grad)=0 and Div(Curl)=0

$$\nabla \times \nabla f=0 ~~ \text{for any} ~~ f$$

$$\nabla \cdot \nabla \times f = 0 ~~ \text{for any} ~~ f$$

## 5. Gauss's Divergence Theorem

The flux of a vector field $\vec F$ through a cloased surface "S" is equal to the integral of the divergence $\vec \nabla \cdot $ over the entire enclosed volumne "V".

$$\iiint\limits_V {\vec \nabla  \cdot \vec FdV = }\iint\limits_S {\vec F \cdot \vec n}dS$$

### Mass conservation in volume

Total mass in $V = \iiint_V \rho dV$, such that

$$\frac{d}{dt} \iiint_V \rho dV = - \iint_S \rho \vec F \cdot \vec n dV = - \iiint_V\vec \nabla  \cdot \rho \vec FdV$$ 

lead to 

$$ \iiint_V \left(\frac{\partial\rho}{\partial t}+ \vec \nabla  \cdot \rho \vec F \right) dV = 0 $$

which is Mass Conservative Equation. 
$$\frac{\partial\rho}{\partial t}+ \vec \nabla  \cdot \rho \vec F = 0$$

## 6. Stokes' Theorem and Green's Theorem

### 6.1 Stokes' Theorem

$$\iint_S \left(\vec \nabla \times \vec F\right)\cdot d \vec A = \int_{\partial S} \vec F \cdot d \vec{s} $$


### 6.2 Green's Theorem 

Stokes on 'flat' $S$,

$$\iint_S \left( \frac{\partial F_2}{\partial x} - \frac{\partial F_1}{\partial y}\right) dx dy = \int_{\partial S} F_1dx + F_2 dy $$

Examples for Stokes' Theorem and Green's Theorem.

## 7. Vector Calculas

Q: Are all vector fields $\vec F$ the gradient of some scalar filed (potential) $f$ ?

A: No. Only special Vector Field's $\vec F = \nabla f$ (potential flows).

Means gradient flows $\nabla f$ are curl free. 

| Flow | Grad flows | Potential flows |
| :-----:| :----: | :----: |
| $\vec F$ | $\vec F = \nabla f$ | $\nabla \cdot \vec F =0 $ & $\nabla \times \vec F =0 $ |
|  | irrotational | incompressible & irrotational|

### 7.1 Gradient flows $\vec F = \nabla f$

Gradient field are conservative: no energy gain/lost.

Moving mass/change around closed loop.

### 7.2 Helmholtz decomposition

Generic flow $\vec F = - \vec \nabla \phi + \vec \nabla \times \vec A$

conservative irrotational incompressible

## 8. Potential Flow

Steady, incompressible, irrotational flow $\vec V = \left[ \begin{array}{} V_1 \\ V_2\end{array} \right]$

Steady: $\frac{\partial \vec V}{\partial t} = 0$

incompressible: $\vec \nabla \cdot \vec V = 0$, $\frac{\partial V_1}{x}+\frac{\partial V_2}{\partial y}=0$

irrotational: $\vec \nabla \times \vec V = 0$, $\frac{\partial V_2}{\partial x}-\frac{\partial V_1}{\partial y} = 0$

These satisfied automatically if $\vec V = \nabla \phi$ for a real-valued potential $\phi$ that satisfies, 
$$\vec \nabla^2 \phi = 0$$ 
Laplace's Equation.

### 8.1 Example potential flow

$\Phi (z) = \phi(x,y)+i \psi(x,y)$

complex potential flow = potential flow + steam flow


### 8.2 Laplace's Equation

- Zero net cirulations acround 'C'
$$\oint\limits_C {\vec V \cdot \vec u} ds = \iint\limits_{{\text{inside }}C} {\vec \nabla  \times \vec V}dA = 0\ $$

- Zero net flux across 'C'
$$\oint\limits_C {\vec V \cdot \vec n} dS = \iint\limits_{{\text{inside }}C} {\vec \nabla  \cdot \vec V}dA = 0$$

- 1. PDE
- 2. Vector field $\vec V$
- 3. ODE for particle in vector field

## 9. Partial Differential Equations (PDE)

Multivariate functions and their partial derivatives

- Wave Equation
$$\bm u_{tt} = c^2 \nabla^2 \bm{u}$$

- Heat Equation
$$\bm u_{t} = \alpha^2 \nabla^2 \bm{u}$$

- Laplace's Equation
$$\nabla^2 \bm{u}=0$$

Linear = Superposition holds

- Nonlinear Burgers Equation 
$$\bm{u}_t + \bm{u} \bm{u}_x = \nu \bm{u}_{xx}$$

### 9.1 Laplace's Equation

- Gravitation (always from mass)
$\vec F = -\nabla V, ~~ V =-\frac{mMG}{r}$

- Electcostatic potential (Coulumb's law similar to gravitation)

- Heat conduction (steady state)
$\frac{\partial T}{\partial t} = \alpha^2 \nabla^2 T=0$

- Incompressible, irrotational flow
$\vec V = \nabla \phi$

### 9.2 Poisson's Equation

$$\nabla^2 \phi =f $$
where $f$ is a forcing function.

## 10. Heat Equation (A parabolic PDE for Energy Conservation)

### 10.1 Deriving the 1D Heat Equation

Considering the temperature distribution $\bm{u}(x,t)$ in a thin metal rod.

The rate of change of heat energy in time = heat flux through boundary to neigh + heat energy generated in $x,t$.

Heat energy 

$$\frac{\partial}{\partial t}\left( c(x)\rho(x)u(x,t)\right) = -\frac{\partial q(x,t)}{\partial x} + Q(x,t)$$
where $q(x,t)$ is the heat flux from left to right, it is expressed
$$q(x,t) = -\Kappa \frac{\partial u}{\partial x}$$

So
$$c(x)\rho(x)\frac{\partial u}{\partial t} = \Kappa \frac{\partial u}{\partial x} +Q(x,t)$$

where $c$, $\rho$, $\Kappa$ constant in the space. We have

$$\frac{\partial u}{\partial t} = \frac{\Kappa}{c\rho} \frac{\partial^2 u}{\partial x^2} + \frac{1}{c\rho}Q(x,t)$$

means

$$\bm u_{t} = \alpha^2 \nabla^2 \bm{u}$$

### 10.2 Solve


## 11. Seperation of Variables ... to solve Laplace's Eqn $\nabla^2u=0$ in 2D

To solve the problem

$\begin{equation} 
    \nabla^2u=0
\end{equation}$

where the boundary conditions:

$\begin{equation} 
    \begin{array}{}
    u(x,0) = 0 \\
    u(x,H) = 0 \\
    u(0,y) = 0 \\
    u(L,y) = f_3(y)
    \end{array}
\end{equation}$

where $H$, $L$ is the computational domain geometry.

For linear partial differential equation (PDE), we have $u(x,y) = F(x) G(y)$, so

$\begin{equation} 
    0 = \nabla^2u = F_{xx}(x)G(y) + F(x)G_{yy}(y)
\end{equation}$

The common requirement of the PDE is

$\begin{equation} 
    \frac{F_{xx}}{F} = - \frac{G_{yy}}{y} = \lambda^{const.}
\end{equation}$

Seperately, for $G(y)$, since there is no boundary condition

$\begin{equation} 
    G_{yy} = -\lambda^{const.} G
\end{equation}$

where $G(0)=G(H)=0$, so, $G(y)=\sin (\sqrt{\lambda} y)$, then $\sqrt{\lambda}H=n\pi$, such that

$\begin{equation} 
    \lambda = \left( \frac{n\pi}{H} \right)^2, ~~n=1,2,3,\dots
\end{equation}$

Finally we get

$\begin{equation} 
    G(y) = \sin \left(\frac{n\pi}{H} y \right), ~~n=1,2,3,\dots
\end{equation}$

Let solve the $F(x)$ with boundaries, we already get some prerequirement of $\lambda$, so

$\begin{equation} 
    F_{xx} = \left(\frac{n\pi}{H}\right)^2 F
\end{equation}$

with the assumption $F(x) = A_n e^{\frac{n\pi}{H}x} + B_n e^{-\frac{n\pi}{H}x}$

When $x=0$, $F(0) = A_n +B_n =0$, leads to $B_n = -A_n$, such that

$\begin{equation} 
    F(x) = A_n e^{\frac{n\pi}{H}x} - A_n e^{-\frac{n\pi}{H}x} = 2A_n \sinh \left({\frac{n\pi}{H}x} \right)
\end{equation}$

So, we have finally solution of Laplace equation (for 1 boundary condition)

$\begin{equation} 
    u(x,y) = F(x)G(y) = \sum_{n=1}^\infty 2A_n\left[ \sinh \left( {\frac{n\pi}{H}x} \right) \right] \left[ \sin \left(\frac{n\pi}{H} y \right) \right]
\end{equation}$

When we get all the solutions (4 boundary condtions for this problem) for the every single boundary condition. Note, here the assumption is $\lambda > 0$.

For another boundary condition, such as
$\begin{equation}
    u(L,y) = f(y) = \sum_{n=1}^\infty A_n\left[ \sinh \left( {\frac{n\pi}{H}L} \right) \right] \left[ \sin \left(\frac{n\pi}{H} y \right) \right]
\end{equation}$

Set one solution,

$\begin{equation}\begin{array}{ll}
    \int_0^H f(y) \sin \left(\frac{m\pi}{H} y \right) dy &= \int_0^HA_m\left[ \sinh \left( {\frac{m\pi}{H}L} \right) \right] \left[ \sin \left(\frac{m\pi}{H} y \right) \sin \left(\frac{m\pi}{H} y \right) \right] \\
    &= \frac{A_mH}{2} \sinh \left( {\frac{m\pi}{H}L} \right)
\end{array}\end{equation}$

So, 
$A_m = \frac{2}{H  \sinh \left( {\frac{m\pi}{H}L} \right)} \int_0^H f(y) \sin \left(\frac{m\pi}{H} y \right) dy$

## 13. The Wave Equation

### 13.1 Deriving the Wave Equation

$\begin{equation} u_{tt} = c^2 u_{xx} \end{equation}$

where $c$ is wave speed, $u(x,t)$ satisfy boundaries, like, $u(0,t) = u(L,t)=0$.

Assumptions

* Graivity is neglegable
* deflections are small
* tension same throughout $x$
* $\rho  \Delta x$ is mass

For each finite length of string,

$\begin{equation} ma = \rho \Delta x u_{tt} \end{equation}$

$\begin{equation} F = T  \sin(\theta+\delta \theta) -T \sin{\theta} = T\left[ \frac{\partial u}{\partial x}(x+\Delta x,t) - \frac{\partial u}{\partial x}(x,t) \right] \end{equation}$

$\begin{equation} u_{tt} = \frac{T}{\rho \Delta x} \left[ \frac{\partial u}{\partial x}(x+\Delta x,t) - \frac{\partial u}{\partial x}(x,t) \right] = \lim_{\Delta x \to 0}  \frac{T}{\rho } \frac{1}{\Delta x} \left[ \frac{\partial u}{\partial x}(x+\Delta x,t) - \frac{\partial u}{\partial x}(x,t) \right] = \frac{T}{\rho } \frac{\partial^2u}{\partial x^2} \end{equation}$

where $c^2 = \frac{T}{\rho }$.

### 13.2 Solving the Wave Equation with Seperation of Variables

* 1. $u(x,t) = F(x)G(t)$
* 2. $u_{tt} = F(x)G^{\prime\prime}(t)$, $u_{xx} = F^{\prime\prime}(x)G(t)$, 
    such that, $F(x)G^{\prime\prime}(t)= c^2 F^{\prime\prime}(x)G(t)$, then, $\begin{equation} \frac{F^{\prime\prime}}{F(x)} = \frac{1}{c^2} \frac{G^{\prime\prime}}{G(t)} = -\lambda^2\end{equation}$

* 3. $\frac{F^{\prime\prime}}{F(x)}=-\lambda^2$ leads to $F^{\prime\prime}+\lambda^2F = 0$, solution is $\begin{equation} F(x) = \beta_n \sin{\frac{n\pi}{L}x} \end{equation}$
* 4. $G^{\prime\prime} = -c^2\lambda^2G$ leads to $G^{\prime\prime} + c^2\lambda^2G = 0$, solution is $\begin{equation} G(t) = k_n \cos{\frac{n\pi c}{L}t} \end{equation}$
* 5. Combine the solutions $\begin{equation} u(x,t) = \sum_{n=0}^{\infty} c_n \sin{\left(\frac{n\pi}{L}x\right)} \cos{\left(\frac{n\pi c}{L}t\right)} \end{equation}$

This is one of Hyperbolic PDE.

As a person steped on a slask line, as a example

$\begin{equation} u_{tt} = c^2 u_{xx} - mg \delta(x-x_0) \end{equation}$

### 13.3 The method of Characteristics and Wave Motion

$\begin{equation}
    u(x,t) = f(x+ct)+f(x-ct)
\end{equation}$

To confirm the $f(x+ct)$ and $f(x-ct)$ are solutions of u(x,t),

$\begin{equation}
    f_{tt} = c^2 f^{\prime\prime}(x+ct)
\end{equation}$

Substitute into Wave Equation (r13), it is confirmed.

