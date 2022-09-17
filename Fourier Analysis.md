<script type="text/javascript"
src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# Fourier Analysis

## 1. Fourier Analysis Overview (Wavelets)

### 1.1 Overview

(Date-driven Science and Engineering Chapter 2)

as kind of
Coordinate transform

$u(x,y,t)$

$u_t=\alpha \nabla^2 u $

SVD = Data-driven FFT

Hilbert Space

Fast Fourier Transform (FFT)

### 1.2 Fourier Series

Define,

$\begin{equation} \left<f(x),g(x) \right> = \int_a^bf(x)\bar g(x)dx \end{equation}$

This is a kind of projection operation takes the function $f(x)$ project into many Sine or Cosine orthogonal coordinates. 

$\begin{equation} \vec f = \left< \vec f, \vec x \right>\frac{\vec x}{\Vert \vec x \Vert^2} + \left< \vec f, \vec y \right>\frac{\vec y}{\Vert \vec y \Vert^2} = \left< \vec f, \vec u \right>\frac{\vec u}{\Vert \vec u \Vert^2} + \left< \vec f, \vec v \right>\frac{\vec v}{\Vert \vec v \Vert^2}\end{equation}$

As we all know,

$\begin{equation} f(x)=\frac{A_0}{2} + \sum_{k=1}^{\inf}\left( A_k \cos (kx) + B_k \sin(kx) \right) \end{equation}$

such that, 

$\begin{equation} A_k=\frac{1}{\pi} \int_{-\pi}^{\pi}f(x)\cos(kx)dx = \frac{1}{\left\Vert cos(kx) \right\Vert^2} \left< f(x),cos(kx)\right>\end{equation}$

$\begin{equation} B_k=\frac{1}{\pi} \int_{-\pi}^{\pi}f(x)\sin(kx)dx = \frac{1}{\left\Vert sin(kx) \right\Vert^2} \left< f(x),sin(kx)\right>\end{equation}$

Found the boundary of the integral taking $[0, L]$, $f(x)$ changes into,

$\begin{equation} f(x)=\frac{A_0}{2} + \sum_{k=1}^{\inf}\left( A_k \cos \left( \frac{2\pi k x}{L} \right) + B_k \sin \left( \frac{2\pi k x}{L} \right) \right) \end{equation}$

Such that, 

$\begin{equation} A_k=\frac{2}{L} \int_{0}^{L}f(x)\cos \left( \frac{2\pi k x}{L} \right)dx \end{equation}$

$\begin{equation} B_k=\frac{2}{L} \int_{0}^{L}f(x)\sin\left( \frac{2\pi k x}{L} \right)dx \end{equation}$

### 1.3 Inner Product of Hilbert Space

#### Inner product

$\begin{equation} \left<f(x),g(x) \right> = \int_a^bf(x)\bar g(x)dx \end{equation}$

For vectors $\underline f = \left[ f_1 \quad f_2 \cdots f_k \cdots f_n \right]^T $ and $\underline g = \left[ g_1 \quad g_2 \cdots g_k \cdots g_n \right]^T$, we have their inner product

$\begin{equation} \left< \underline f, \underline g \right> =\underline g ^T \underline f = \sum_{k=1}^{n} f_k g_k \end{equation}$

#### Riemann approximation

$\begin{equation} \left< \underline f, \underline g \right>\Delta x =\sum_{k=1}^{n} f(x_k) g(x_k)\Delta x \end{equation}$

where $\Delta x = L/(n-1) = (b-a)/(n-1)$.

### 1.4 Complex Fourier Series

