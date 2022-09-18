<script type="text/javascript"
src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

## Vector Calculus and Partial Differential Equations

### Big Picture Overview

Conservation Laws
* Mass
* Momentum
* Energy

Vector Calc
* $\nabla $   grad
* $\nabla \cdot$  div
* $\nabla \times$ curl


### Seperation of Variables ... to solve Laplace's Eqn $\nabla^2u=0$ in 2D

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

## The Wave Equation

### Deriving the Wave Equation

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

### Solving the Wave Equation with Seperation of Variables

* 1. $u(x,t) = F(x)G(t)$
* 2. $u_{tt} = F(x)G^{\prime\prime}(t)$, $u_{xx} = F^{\prime\prime}(x)G(t)$, 
    such that, $F(x)G^{\prime\prime}(t)= c^2 F^{\prime\prime}(x)G(t)$, then, $\begin{equation} \frac{F^{\prime\prime}}{F(x)} = \frac{1}{c^2} \frac{G^{\prime\prime}}{G(t)} = -\lambda^2\end{equation}$

* 3. $\frac{F^{\prime\prime}}{F(x)}=-\lambda^2$ leads to $F^{\prime\prime}+\lambda^2F = 0$, solution is $\begin{equation} F(x) = \beta_n \sin{\frac{n\pi}{L}x} \end{equation}$
* 4. $G^{\prime\prime} = -c^2\lambda^2G$ leads to $G^{\prime\prime} + c^2\lambda^2G = 0$, solution is $\begin{equation} G(t) = k_n \cos{\frac{n\pi c}{L}t} \end{equation}$
* 5. Combine the solutions $\begin{equation} u(x,t) = \sum_{n=0}^{\infty} c_n \sin{\left(\frac{n\pi}{L}x\right)} \cos{\left(\frac{n\pi c}{L}t\right)} \end{equation}$

This is one of Hyperbolic PDE.

As a person steped on a slask line, as a example

$\begin{equation} u_{tt} = c^2 u_{xx} - mg \delta(x-x_0) \end{equation}$

### The method of Characteristics and Wave Motion

$\begin{equation}
    u(x,t) = f(x+ct)+f(x-ct)
\end{equation}$

To confirm the $f(x+ct)$ and $f(x-ct)$ are solutions of u(x,t),

$\begin{equation}
    f_{tt} = c^2 f^{\prime\prime}(x+ct)
\end{equation}$

Substitute into Wave Equation (13), it is confirmed.

