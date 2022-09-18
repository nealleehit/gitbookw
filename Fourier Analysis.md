@[TOC](目录)

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

$\begin{equation} \left<f(x),g(x) \right> = \int_{-\pi}^{\pi} f(x)\bar g(x)dx \end{equation}$

where $f(x) = \sum_{k=-\infty}^{\infty} c_k e^{ikx}$ and $e^{ikx} = \cos(kx) +i \sin(kx)$.

$\begin{equation} f(x) = \sum_{k=-\infty}^{\infty} c_k e^{ikx} = \sum_{k=-\infty}^{\infty} (\alpha_k+i \beta_k)\left[\cos(kx) +i \sin(kx) \right]\end{equation}$

Define $\psi_k :=e^{ikx}$ as the orthogonal base vectors. So we have

$\begin{equation} \left< \psi_j, \psi_k \right> 
= \int_{-\pi}^{\pi}e^{ijx} e^{-ikx}dx 
=\int_{-\pi}^{\pi} e^{i(j-k)x}dx 
= \frac{1}{i(j-k)}\left[ e^{i(j-k)x}\right]^\pi_{-\pi} 
= \left\{ \begin{array}{lc} 0 ~~\rm{ if }~~ j \ne k \\ 2\pi ~~\rm{if } ~~j= k \end{array} \right.\end{equation}$

Then $f(x)$ can be writen as

$$f(x) = \frac{1}{2\pi} \sum_{k=-\infty}^{\infty} \underbrace {\left< {f\left( x \right),{\psi _k}} \right> }_{c_k}\underbrace {\psi _k}_{e^{ikx}}$$

Ex. Inner Product

```Matlab
clear all, close all, clc

f = [0 0 .1 .2  .25  .2 .25 .3 .35 .43 .45 .5 .55  .5 .4 .425 .45 .425 .4 .35 .3 .25 .225 .2 .1 0 0];
g = [0 0 .025 .1  .2  .175 .2 .25 .25 .3 .32 .35 .375  .325 .3 .275 .275 .25 .225 .225 .2 .175 .15 .15 .05 0 0] -0.025;

x = 0.1*(1:length(f));

xf = (.01:.01:x(end));
ff = interp1(x,f,xf,'cubic')

gf = interp1(x,g,xf,'cubic')


plot(xf(20:end-10),ff(20:end-10),'k','LineWidth',1.5)
hold on
plot(x(2:end-1),f(2:end-1),'bo','MarkerFace','b')
plot(xf(20:end-10),gf(20:end-10),'k','LineWidth',1.5)
plot(x(2:end-1),g(2:end-1),'ro','MarkerFace','r')


xlim([.1 2.7])
ylim([-.1 .6])
set(gca,'XTick',[.2:.1:2.6],'XTickLabels',{},'LineWidth',1.2)
set(gca,'YTick',[]);
box off

set(gcf,'Position',[100 100 550 250])

set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/InnerProduct');

% %%
% xc = x;
% fc = f;
% n = length(x);
% hold on
% fapx = 0*ff;
% dx = xc(2)-xc(1);
% L = xc(end)-xc(1);
% L = 2.5
% A0 = (1/pi)*sum(fc.*ones(size(xc)))*dx*L;
% fapx = fapx + A0/2;
% for k=1:10
%     Ak = (1/pi)*sum(fc.*cos(2*pi*k*xc/L))*dx*L;
%     Bk = (1/pi)*sum(fc.*sin(2*pi*k*xc/L))*dx*L;
% 
%     fapx = fapx + Ak*cos(2*k*pi*xf/L) + Bk*sin(2*k*pi*xf/L);
% end
%     plot(xf,fapx,'k')
```

Ex2. Fourier Sines Products

```Matlab
clear all, close all, clc

kmax = 7;

dx = 0.001;
L = pi;
x = (-1+dx:dx:1)*L;
f = 0*x;
n = length(f);
nquart = floor(n/4);
nhalf = floor(n/2);

f(nquart:nhalf) = 4*(1:nquart+1)/n;
f(nhalf+1:3*nquart) = 1-4*(0:nquart-1)/n;
subplot(3,1,1)
plot(x,f,'-','Color',[0 0 0],'LineWidth',1.5)
ylim([-.2 1.5])
xlim([-1.25*L 1.25*L])
set(gca,'LineWidth',1.2)
set(gca,'XTick',[-L 0 L],'XTickLabels',{});%{'-L','0','L','2L'})
set(gca,'YTick',[0 1],'YTickLabels',{});
box off

CC = colormap(jet(8));
% CCsparse = CC(5:5:end,:);
% CCsparse(end+1,:) = CCsparse(1,:);
CCsparse = CC(1:3:end,:);
%
subplot(3,1,2)
L = pi;
A0 = sum(f.*ones(size(x)))*dx;
plot(x,A0+0*f,'-','Color',CC(1,:)*.8,'LineWidth',1.2);
hold on
fFS = A0/2;
for k=1:kmax
    A(k) = sum(f.*cos(pi*k*x/L))*dx;
    B(k) = sum(f.*sin(pi*k*x/L))*dx;
    plot(x,A(k)*cos(k*pi*x/L),'-','Color',CC(k,:)*.8,'LineWidth',1.2);
%     plot(x,B(k)*sin(2*k*pi*x/L),'k-','LineWidth',1.2);
    fFS = fFS + A(k)*cos(k*pi*x/L) + 0*B(k)*sin(k*pi*x/L);
end
ylim([-.7 .7])
xlim([-1.25*L 1.25*L])
set(gca,'LineWidth',1.2)
set(gca,'XTick',[-L 0 L],'XTickLabels',{});%{'-L','0','L','2L'})
set(gca,'YTick',[-.5 0 .5],'YTickLabels',{});
box off
% 
subplot(3,1,1)
hold on
plot(x,fFS,'-','Color',CC(7,:)*.8,'LineWidth',1.2)
l1=legend('     ','    ')
set(l1,'box','off');
l1.FontSize = 16;


subplot(3,1,3)
A0 = sum(f.*ones(size(x)))*dx;
plot(x,A0+0*f,'-','Color',CC(1,:),'LineWidth',1.2);
hold on
fFS = A0/2;
for k=1:7
    Ak = sum(f.*cos(pi*k*x/L))*dx;
    Bk = sum(f.*sin(pi*k*x/L))*dx;
    plot(x,Ak*cos(k*pi*x/L),'-','Color',CC(k,:)*.8,'LineWidth',1.2);
%     plot(x,Bk*sin(2*k*pi*x/L),'k-','LineWidth',1.2);
    fFS = fFS + Ak*cos(k*pi*x/L) + 0*Bk*sin(k*pi*x/L);
end
ylim([-.06 .06])
xlim([-1.25*L 1.25*L])
set(gca,'LineWidth',1.2)
set(gca,'XTick',[-L 0 L],'XTickLabels',{});%{'-L','0','L','2L'})
set(gca,'YTick',[-.05 0 .05],'YTickLabels',{});
box off

set(gcf,'Position',[100 100 550 400])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/FourierTransformSines');

%% Plot amplitudes
clear ERR
clear A
fFS = A0/2;
A(1) = A0/2;
ERR(1) = norm(f-fFS);
kmax = 100;
for k=1:kmax
    A(k+1) = sum(f.*cos(2*pi*k*x/L))*dx*2/L;
    B(k+1) = sum(f.*sin(2*pi*k*x/L))*dx*2/L;
%     plot(x,B(k)*sin(2*k*pi*x/L),'k-','LineWidth',1.2);
    fFS = fFS + A(k+1)*cos(2*k*pi*x/L) + 0*B(k+1)*sin(2*k*pi*x/L);
    ERR(k+1) = norm(f-fFS)/norm(f);
end
thresh = median(ERR)*sqrt(kmax)*4/sqrt(3);
r = max(find(ERR>thresh));
r = 7;
subplot(2,1,1)
semilogy(0:1:kmax,A,'k','LineWidth',1.5)
hold on
semilogy(r,A(r+1),'bo','LineWidth',1.5)
xlim([0 kmax])
ylim([10^(-7) 1])
subplot(2,1,2)
semilogy(0:1:kmax,ERR,'k','LineWidth',1.5)
hold on
semilogy(r,ERR(r+1),'bo','LineWidth',1.5)
xlim([0 kmax])
ylim([3*10^(-4) 20])
set(gcf,'Position',[100 100 500 300])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', '../figures/FourierTransformSinesERROR');
```

Ex 3. Fourier Sines

```Matlab
clear all, close all, clc

% Define domain
dx = 0.001;
L = pi;
x = (-1+dx:dx:1)*L;
n = length(x);   nquart = floor(n/4);

% Define hat function
f = 0*x;
f(nquart:2*nquart) = 4*(1:nquart+1)/n;
f(2*nquart+1:3*nquart) = 1-4*(0:nquart-1)/n;
plot(x,f,'-k','LineWidth',1.5), hold on

% Compute Fourier series
CC = jet(20);
A0 = sum(f.*ones(size(x)))*dx;
fFS = A0/2;
for k=1:20
    A(k) = sum(f.*cos(pi*k*x/L))*dx; % Inner product
    B(k) = sum(f.*sin(pi*k*x/L))*dx;
    fFS = fFS + A(k)*cos(k*pi*x/L) + B(k)*sin(k*pi*x/L);
    plot(x,fFS,'-','Color',CC(k,:),'LineWidth',1.2)
end


%% Plot amplitudes
figure
clear ERR
clear A
fFS = A0/2;
A(1) = A0/2;
ERR(1) = norm(f-fFS);
kmax = 100;
for k=1:kmax
    A(k+1) = sum(f.*cos(pi*k*x/L))*dx;
    B(k+1) = sum(f.*sin(pi*k*x/L))*dx;
%     plot(x,B(k)*sin(2*k*pi*x/L),'k-','LineWidth',1.2);
    fFS = fFS + A(k+1)*cos(k*pi*x/L) + B(k+1)*sin(k*pi*x/L);
    ERR(k+1) = norm(f-fFS)/norm(f);
end
thresh = median(ERR)*sqrt(kmax)*4/sqrt(3);
r = max(find(ERR>thresh));
r = 7;
subplot(2,1,1)
semilogy(0:1:kmax,A,'k','LineWidth',1.5)
hold on
semilogy(r,A(r+1),'bo','LineWidth',1.5)
xlim([0 kmax])
ylim([10^(-7) 1])
subplot(2,1,2)
semilogy(0:1:kmax,ERR,'k','LineWidth',1.5)
hold on
semilogy(r,ERR(r+1),'bo','LineWidth',1.5)

```

### 1.5 Gibbs Phenomena

Ex Gibbs

```Matlab
clear all, close all, clc
dx = 0.01;  L = 10;
x = 0:dx:L;
n = length(x); nquart = floor(n/4);

f = zeros(size(x));
f(nquart:3*nquart) = 1;

A0 = sum(f.*ones(size(x)))*dx*2/L;
fFS = A0/2;
for k=1:100
    Ak = sum(f.*cos(2*pi*k*x/L))*dx*2/L;
    Bk = sum(f.*sin(2*pi*k*x/L))*dx*2/L;
    fFS = fFS + Ak*cos(2*k*pi*x/L) + Bk*sin(2*k*pi*x/L);
end   

plot(x,f,'k','LineWidth',2), hold on
plot(x,fFS,'r-','LineWidth',1.2)
```

Ex Gibbs_production
```Matlab
clear all, close all, clc

dx = 0.01;
L = 10;
x = 0:dx:L;

f = zeros(size(x)); 
f(floor(length(f)/4):floor(3*length(f)/4)) = 1;%+f(floor(length(f)/4):floor(3*length(f)/4));
% fFS = zeros(size(x));

A0 = sum(f.*ones(size(x)))*dx*2/L;
for m=100%1:100
    fFS = A0/2;
    for k=1:m
        Ak = sum(f.*cos(2*pi*k*x/L))*dx*2/L;
        Bk = sum(f.*sin(2*pi*k*x/L))*dx*2/L;
        fFS = fFS + Ak*cos(2*k*pi*x/L) + Bk*sin(2*k*pi*x/L);
    end
    
    plot(x,f,'k','LineWidth',2)
    hold on
    plot(x,fFS,'r-','LineWidth',1.2)
    
    pause(0.1)
end

set(gca,'XTick',[],'XTickLabels',{},'LineWidth',1.2)
set(gca,'YTick',[]);
box off
% axis off

set(gcf,'Position',[100 100 550 200])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', '../figures/Gibbs');
```

Ex Gibbs_Movie
```Matlab
clear all, close all, clc

dx = 0.01;
L = 2*pi;
x = 0:dx:L;

f = zeros(size(x));% sin(x/L);%ones(size(x));
f(floor(length(f)/4):floor(3*length(f)/4)) = 1+f(floor(length(f)/4):floor(3*length(f)/4));
% f = (x-pi).^2;
fFS = zeros(size(x));

A0 = (1/pi)*sum(f.*ones(size(x)))*dx;
for m=1:100
    fFS = A0/2;
    for k=1:m
        Ak = (1/pi)*sum(f.*cos(2*pi*k*x/L))*dx;
        Bk = (1/pi)*sum(f.*sin(2*pi*k*x/L))*dx;
        fFS = fFS + Ak*cos(2*k*pi*x/L) + Bk*sin(2*k*pi*x/L);
    end
    
    plot(x,f,'k')
    hold on
    plot(x,fFS,'r-')
    
    pause(0.1)
end
```

## 2. Fourier Transform Topics

### 2.1 Fourier Transforms 

If we already have the complex Fourier expression,

$\begin{equation} f(x) = \sum_{k=-\infty}^{\infty} c_k e^{ik\pi x/L} \end{equation}$

where $\omega_k = \frac{k\pi}{L} = k \Delta \omega$, $\Delta \omega = \pi/L$, and $c_k = \frac{1}{2\pi} \left< f(x), \psi_k \right>  =\frac{1}{2L}\int_{-L}^{L} f(x) \underbrace{e^{-ik\pi x/L}}_{\psi_k}  dx $. The Fourier Transform expressed as,

$\begin{equation} f(x) = \lim_{\Delta \omega \to 0}\sum_{k=-\infty}^{\infty} \frac{\Delta \omega}{2\pi} \int_{-\pi/\Delta \omega}^{\pi/\Delta \omega} f(\xi) e^{-ik \Delta \omega \xi} d\xi e^{ik \Delta \omega x} \\ = \int_{-\infty}^{\infty} \frac{1}{2\pi} \int_{-\infty}^{\infty} f(\xi) e^{-i \omega \xi} d\xi e^{i \omega x} d\omega\end{equation}$

So the Fourier Transforms,

$\begin{equation} \hat f(\omega) = \cal{F}\left(f(x)\right) = \int_{-\infty}^{\infty} f(x) e^{-i \omega x} dx \end{equation}$

$\begin{equation} f(x) = \cal{F}^{-1} \left( \hat f (\omega)\right) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat f(\omega) e^{i \omega x} d\omega \end{equation}$

These are the Fourier Transform Pair.

### 2.2 Fourier Transform and Derivatives

$\begin{equation} 
    \begin{array}{ll} 
    \cal{F}\left( \frac{d}{dx}f(x)\right) &= \int_{-\infty}^{\infty} \underbrace{\frac{df(x)}{dx}}_{dv} \underbrace{e^{-i \omega x}}_{u} dx \\ 
    & = \underbrace{\left[ f(x)e^{-i \omega x} \right]^{\infty}_{-\infty}}_{uv=0} - \int_{-\infty}^{\infty} \underbrace{df(x)}_{v} \underbrace{(-i\omega e^{-i \omega x})}_{du} dx \\ 
    & = i \omega \underbrace{\int_{-\infty}^{\infty} f(x) e^{-i \omega x} dx}_{\cal{F}\left( f(x)\right)} \\ 
    & = i \omega \cal{F}\left( f(x)\right) 
    \end{array} 
\end{equation}$

### 2.3 Fourier Transform and Convolution Integral

Convolution integral

$\begin{equation} 
    \left( f * g \right) = \int_{-\infty}^{\infty} f(x-\xi) g(\xi) d\xi 
\end{equation}$

Fourier 

$\begin{equation} 
    \cal{F}\left( f*g \right) =  \cal{F}\left( f \right) \cal{F}\left( g \right) = \hat f \hat g 
\end{equation}$

While

$\begin{equation} 
    \begin{array}{ll}
    \cal{F}^{-1}\left( \hat f \hat g \right)(x) & = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat f(\omega) \hat g(\omega) e^{i\omega x} d\omega \\
    & = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat f(\omega)\left( \int_{-\infty}^{\infty} g(y) e^{-i\omega y} dy \right)e^{i\omega x} d\omega \\
    & = \frac{1}{2\pi} \int_{-\infty}^{\infty} g(y) \int_{-\infty}^{\infty} \hat f(\omega)  e^{i\omega (x- y)} d\omega dy \\
    & = \int_{-\infty}^{\infty} g(y) f(x-y) dy \\
    & = f * g
    \end{array}
\end{equation}$

### 2.4 Parseval's Theorem

$$ \cal{F} \left( \alpha f(x) +\beta g(x)\right) = \alpha \cal{F}(f) + \beta \cal{F}(g)$$

The Parseval's Theorem is 

$\begin{equation}
    \int_{-\infty}^{\infty} \left\vert \hat f(\omega) \right\vert^2 d\omega = 2\pi\int_{-\infty}^{\infty} \left\vert f(x) \right\vert^2 dx
\end{equation}$

## 3. The Discrete Fourier Transform (DFT)

### 3.1 DFT

$\begin{equation} \begin{array}{}
    \hat f_k & = \sum_{j=0}^{n-1} \left( f_j e^{-i2\pi jk/n} \right)\\
    f_k & = \frac{1}{n} \sum_{j=0}^{n-1} \left( \hat f_j e^{i2\pi jk/n} \right)
\end{array}\end{equation}$

such that, $\omega_n = e^{-2\pi i/n}$

$\begin{equation} 
    \left[ \begin{array}{}
    \hat f_0 \\ \hat f_1 \\ \hat f_2 \\ \vdots \\ f_n
    \end{array} \right] = 
    \left[ \begin{array}{c} 
    1 & 1 & 1 & \cdots & 1 \\
    1 & \omega_n & \omega_n^2 & \cdots & \omega_n^{n-1} \\
    1 & \omega_n^2 & \omega_n^4 & \cdots & \omega_n^{2(n-1)} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    1 & \omega_n^{n-1} & \omega_n^{2(n-1)} & \cdots & \omega_n^{(n-1)^2}
    \end{array} \right] 
    \left[ \begin{array}{} f_0 \\ f_1 \\ f_2 \\ \vdots \\ f_n \end{array} \right] 
\end{equation}$
where $\omega_n = e^{-2\pi i/n}$. This matrix expression is the DFT matrix.

```Matlab
clear all, close all, clc
n = 256;
w = exp(-i*2*pi/n);

% Slow
for i=1:n
    for j=1:n
        DFT(i,j) = w^((i-1)*(j-1));
    end
end

% Fast
[I,J] = meshgrid(1:n,1:n);
DFT = w.^((I-1).*(J-1));
imagesc(real(DFT))
```

## 4. Fast Fourier Transform (FFT)

The computational cost for DFT is $\cal{O}(n^2)$, for FFT is$\cal{O}(n\log (n))$, Which means for a series of data $n=4.4 \times 10^5$, the computation for DFT and FFT is $10^{11}$ and $10^6$, respectively.

### 4.1 The FFT Algorithm

As in Eqn. $(25)$, 

$\begin{equation} 
    \hat{\mathbf{f}}  = \cal{F}_{1024} \mathbf{f} = \left[ \begin{array}{}
    I_{512} & -D_{512} \\ I_{512} & -D_{512}\end{array} \right] 
    \left[ \begin{array}{c} 
    \cal{F}_{512} & \cal{0} \\
    \cal{0}  & \cal{F}_{512} 
    \end{array} \right] \left[ \begin{array}{}
    \mathbf{f}_{even} \\ \mathbf{f}_{odd} \end{array} \right]
\end{equation}$
where $I_{512}$ is identity matrix with dimensionof $512 \times 512$, and
$\begin{equation}
    D_{512} = \left[ \begin{array}{c} 
    1 & 0  &  0 & \cdots & 0 \\
    0 & \omega_n & 0 & \cdots & 0 \\
    0 & 0& \omega_n^2 & \cdots & 0 \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & 0 & \cdots & \omega_n^{511}
    \end{array} \right] 
\end{equation}$

So we can 

$\cal{F}_{1024} \to \cal{F}_{512} \to \cal{F}_{256} \to \cal{F}_{128} \cdots \to \cal{F}_{4} \to \cal{F}_2$.

- Ex. Denoise Data with FFT

```Matlab
clear all

%% create a simple signal with two frequencies
dt = .001;
t = 0:dt:1;
x = sin(2*pi*50*t) + sin(2*pi*120*t);
y = x + 2.5*randn(size(t));  %  add some noise


%% Compute the Fast Fourier Transform FFT
N = length(t);
Y = fft(y,N);  % computes the (fast) discrete fourier transform
PSD = Y.*conj(Y)/N;  % Power spectrum (how much power in each freq)
freq = 1/(dt*N)*(0:N);  %create the x-axis of frequencies in Hz
L = 1:floor(N/2);  % only plot the first half of freqs


%% Use the PSD to filter out noise
indices = PSD>100;   % Find all freqs with large power
PSDclean = PSD.*indices;  % Zero out all others

Y = indices.*Y;  % zero out small Fourier coefficients in Y
yfilt = ifft(Y);     % inverse FFT to get filtered time-domain signal


%% PLOTS
subplot(3,1,1)
plot(-1000,-1000,'k','LineWidth',1.5)
hold on
plot(-1000,-1000,'r','LineWidth',1.2)
plot(t,y,'r','LineWidth',1.2)
plot(t,x,'k','LineWidth',1.5)
axis([0 .25 -5 5])
legend('Clean','Noisy')
set(gca,'LineWidth',1.2,'FontSize',12)

subplot(3,1,3)
plot(t,x,'k','LineWidth',1.5)
hold on
plot(t,yfilt,'b','LineWidth',1.2)
axis([0 .25 -5 5])
legend('Clean','Filtered')
set(gca,'LineWidth',1.2,'FontSize',12)

subplot(3,1,2)
plot(freq(L),PSD(L),'r','LineWidth',1.5)
% xlabel('Frequency (Hz)')
hold on
plot(freq(L),PSDclean(L),'-b','LineWidth',1.2)
legend('Noisy','Filtered')
set(gca,'LineWidth',1.2,'FontSize',12)

set(gcf,'Position',[100 100 550 450])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/DENOISE');
```

* Ex. Computing Derivatives with FFT

```Matlab
clear all, close all, clc

n = 128;
L = 30;
dx = L/(n);
x = -L/2:dx:L/2-dx;
f = cos(x).*exp(-x.^2/25);                    % Function
df = -(sin(x).*exp(-x.^2/25) + (2/25)*x.*f);  % Derivative

%% Approximate derivative using finite Difference...
for kappa=1:length(df)-1
    dfFD(kappa) = (f(kappa+1)-f(kappa))/dx;
end
dfFD(end+1) = dfFD(end);

%% Derivative using FFT (spectral derivative)
fhat = fft(f);
kappa = (2*pi/L)*[-n/2:n/2-1];
kappa = fftshift(kappa);  % Re-order fft frequencies
dfhat = i*kappa.*fhat;
dfFFT = real(ifft(dfhat));

%% Plotting commands
plot(x,df,'k','LineWidth',1.5), hold on
plot(x,dfFD,'b--','LineWidth',1.2)
plot(x,dfFFT,'r--','LineWidth',1.2)
legend('True Derivative','Finite Diff.','FFT Derivative')
```

## 5. Solve PDEs with the FFT

* Ex1. FFT solve Heat PDE

```Matlab
clear all, close all, clc

% Define spatial domain
c = 2;              % Wave speed
L = 20;             % Length of domain 
N = 1000;           % Number of discretization points
dx = L/N;
x = -L/2:dx:L/2-dx; % Define x domain

% Define discrete wavenumbers
kappa = (2*pi/L)*[-N/2:N/2-1];
kappa = fftshift(kappa');    % Re-order fft wavenumbers

% Initial condition 
u0 = sech(x);        
uhat0 = fft(u0);

% Simulate in Fourier frequency domain
dt = 0.025;
t = 0:dt:100*dt;
[t,uhat] = ode45(@(t,uhat)rhsWave(t,uhat,kappa,c),t,uhat0);

% Alternatively, simulate in spatial domain
[t,u] = ode45(@(t,u)rhsWaveSpatial(t,u,kappa,c),t,u0);

% Inverse FFT to bring back to spatial domain
for k = 1:length(t)
    u(k,:) = ifft(uhat(k,:));
end

% Plot solution in time
subplot(1,2,1)
h=waterfall(real(u(1:10:end,:)));
set(h,'LineWidth',2,'FaceAlpha',.5);
colormap(jet/1.5)

subplot(1,2,2)
imagesc(flipud(real(u)));
colormap jet

%% FIGURES (PRODUCTION)
figure
CC = colormap(jet(101));
dt = 0.025;
for k = 1:length(t)
    u(k,:) = real(ifft(uhat(k,:)));
    if(mod(k-1,10)==0)
        plot(x,u(k,:),'Color',CC(k,:),'LineWidth',1.5)
        hold on, grid on, drawnow
    end   
end
% xlabel('Spatial variable, x')
% ylabel('Temperature, u(x,t)')
axis([-L/2 L/2 -.1 1.1])
set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'Position',[100 100 550 220]);
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../../figures/FFTWave1');

%
figure
subplot(1,2,1)
h=waterfall(u(1:10:end,:));
set(h,'LineWidth',2,'FaceAlpha',.5);
colormap(jet/1.5)
% view(22,29)
view(3,21)
set(gca,'LineWidth',1.5)
set(gca,'XTick',[0 500 1000],'XTickLabels',{})
set(gca,'ZTick',[0 .5 1],'ZTickLabels',{})
set(gca,'YTick',[0 5 10],'YTickLabels',{})

subplot(1,2,2)
imagesc(flipud(u));
set(gca,'LineWidth',1.5)
set(gca,'XTick',[0 500 1000],'XTickLabels',{})
set(gca,'YTick',[0 50 100],'YTickLabels',{})

colormap jet
set(gcf,'Position',[100 100 600 250])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../../figures/FFTWave2');
```

* Ex2. FFT solve Berger PDE

```Matlab
clear all, close all, clc
nu=0.001;   % Diffusion constant

% Define spatial domain
L = 20;             % Length of domain 
N = 1000;           % Number of discretization points
dx = L/N;
x = -L/2:dx:L/2-dx; % Define x domain

% Define discrete wavenumbers
kappa = (2*pi/L)*[-N/2:N/2-1];
kappa = fftshift(kappa');    % Re-order fft wavenumbers

% Initial condition 
u0 = sech(x);        

% Simulate PDE in spatial domain
dt = 0.025;
t = 0:dt:100*dt;
[t,u] = ode45(@(t,u)rhsBurgers(t,u,kappa,nu),t,u0);

% Plot solution in time
subplot(1,2,1)
h=waterfall(real(u(1:10:end,:)));
set(h,'LineWidth',2,'FaceAlpha',.5);
colormap(jet/1.5)
view(5,55)

subplot(1,2,2)
imagesc(flipud(real(u)));
colormap jet

%% FIGURES (PRODUCTION)
figure
CC = colormap(jet(100));
dt = 0.1;
for k = 1:100
    if(mod(k-1,10)==0)
        plot(x,real(u(k,:)),'k','LineWidth',1.5)
%         plot(x,real(u(k,:)),'Color',CC(k,:),'LineWidth',1.5)
        hold on, grid on, drawnow
    end   
end
% xlabel('Spatial variable, x')
% ylabel('Temperature, u(x,t)')
axis([-L/2 L/2 -.1 1.1])
set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'Position',[100 100 550 220]);
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../../figures/FFTBurgers1');

%
figure
subplot(1,2,1)
h=waterfall(real(u(1:10:end,:)));
set(h,'LineWidth',2,'FaceAlpha',.5);
colormap(jet/1.5)
view(5,55)
set(gca,'LineWidth',1.5)
set(gca,'XTick',[0 500 1000],'XTickLabels',{})
set(gca,'ZTick',[0 .5 1],'ZTickLabels',{})
set(gca,'YTick',[0 5 10],'YTickLabels',{})
set(gca,'YLim',[1 11])
set(gca,'ZLim',[-.1 1.1])

subplot(1,2,2)
imagesc(flipud(real(u)));
set(gca,'LineWidth',1.5)
set(gca,'XTick',[0 500 1000],'XTickLabels',{})
set(gca,'YTick',[0 50 100],'YTickLabels',{})

colormap jet
set(gcf,'Position',[100 100 600 250])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../../figures/FFTBurgers2');
```

## 6. The Spectrogram and the Gabor Transform

### 6.1 Gabor Transform

$$G(f) = \hat f_g(t, \omega) = \int_{-\infty}^{\infty}f(\tau)e^{-i\omega \tau}g(t-\tau) d\tau \tag{28} $$

Ex. Spectrogram

```Matlab
clear all, close all, clc

t = 0:0.001:2;
f0 = 50;
f1 = 250;
t1 = 2;
x = chirp(t,f0,t1,f1,'quadratic');
x = cos(2*pi*t.*(f0 + (f1-f0)*t.^2/(3*t1^2)));
% There is a typo in Matlab documentation... 
% ... divide by 3 so derivative amplitude matches frequency 

spectrogram(x,128,120,128,1e3,'yaxis')
colormap jet

set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'Position',[100 100 550 200]);
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../../figures/SPECTROGRAM_CHIRP');
```

### 6.2 Uncertainty Principle

$$\left( \int_{\infty}^{\infty} x^2 \left\vert f(x) \right\vert^2 dx \right)\left( \int_{\infty}^{\infty} \omega^2 \left\vert \hat f(\omega) \right\vert^2 d\omega \right) \ge \frac{1}{16\pi^2} \tag{29}$$

## 7. Wavelet and Multiresolution Analysis

### 7.1 Wavelet

Mother wavelet $\Psi(t)$

$$\Psi_{a,b}(t) = \frac{1}{\sqrt a} \phi\left(\frac{t-b}{a}\right) \tag{30}$$

$$\cal{W}_\Psi(f)(a,b) = \left< f(t),\Psi_{a,b}(t) \right>$$

Haar Wavelet 1910 (1-0 wavelet); Doubechies wavelet; Maxican Hat wavelet; Coiflet wavelet, et al.

- Ex. Haar wavelet

```Matlab
clear all, close all, clc

x = 0:.001:1;

n = length(x);
n2 = floor(n/2);
n4 = floor(n/4);

f10 = 0*x;
f10(1:n2-1) = 1;
f10(n2:end) = -1;

f21 = 0*x;
f21(1:n4-1) = 1;
f21(n4:n2-1) = -1;
f21 = f21*sqrt(2);

f22 = 0*x;
f22(n2:n2+n4-1) = 1;
f22(n2+n4:end) = -1;
f22 = f22*sqrt(2);

x = [-1 0 x 1 2];
f10 = [0 0 f10 0 0];
f21 = [0 0 f21 0 0];
f22 = [0 0 f22 0 0];

subplot(3,1,1)
plot(x,f10,'k','LineWidth',2)
xlim([-.2 1.2])
ylim([-1.2 1.2])
set(gca,'XTick',[0 .25 .5 .75 1])
set(gca,'LineWidth',1.2,'FontSize',12);
subplot(3,1,2)
plot(x,f21,'k','LineWidth',2)
xlim([-.2 1.2])
ylim([-1.75 1.75])
set(gca,'XTick',[0 .25 .5 .75 1])
set(gca,'LineWidth',1.2,'FontSize',12);
subplot(3,1,3)
plot(x,f22,'k','LineWidth',2)
xlim([-.2 1.2])
ylim([-1.75 1.75])
set(gca,'XTick',[0 .25 .5 .75 1])
set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'Position',[100 100 550 350]);
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', '../figures/HAAR');

%%
figure
x = -5:.001:5;

fMexHat = (1-x.^2).*exp(-x.^2/2);
plot(x,fMexHat,'k','LineWidth',2)
```

- Ex. Compression image with FFT

```Matlab
clear all, close all, clc
A = imread('../../CH01_SVD/DATA/dog.jpg');
B = rgb2gray(A);

%% FFT Compression
Bt=fft2(B);    % B is grayscale image from above
Btsort = sort(abs(Bt(:)));  % Sort by magnitude

% Zero out all small coefficients and inverse transform
for keep=[.1 .05 .01 .002];
    thresh = Btsort(floor((1-keep)*length(Btsort)));
    ind = abs(Bt)>thresh;      % Find small indices
    Atlow = Bt.*ind;           % Threshold small indices
    Alow=uint8(ifft2(Atlow));  % Compressed image
    figure, imshow(Alow)      % Plot Reconstruction
end
set(f1,'Position',[100 100 600 800])
set(f2,'Position',[100 100 600 800])
```

- Ex. Compression image with wavelet

```Matlab
clear all, close all, clc
A = imread('../../CH01_SVD/DATA/dog.jpg');
B = rgb2gray(A);

%% Wavelet Compression
[C,S] = wavedec2(B,4,'db1');
Csort = sort(abs(C(:))); % Sort by magnitude

for keep =  [.1 .05 .01 .005]
    thresh = Csort(floor((1-keep)*length(Csort)));
    ind = abs(C)>thresh;
    Cfilt = C.*ind;      % Threshold small indices
    
    % Plot Reconstruction
    Arecon=uint8(waverec2(Cfilt,S,'db1'));
    figure, imagesc(uint8(Arecon))
end
```


## 8. The Laplace Transform

PDE $\to$ ODE

ODE $\to$ algebraic

Control Theory

### 8.1 Derive the Laplace Transform

- For some badly behaved function, such as, $f(t) = e^{\lambda t}$, Fourier Transform is not good choice.

- Solution, multiply $f(t)$ by $e^{-\gamma t}$, so that
$$ f(t)e^{-\gamma t} \to 0 ~~~ {as} ~~~ t \to \infty \tag{31}$$

We have $F(t) = f(t)e^{-\gamma t}$, and we can constructe Heaviside function, $H(t) = \left\{ \begin{array}{}0 & t<0 \\ f(t)e^{-\gamma t} & t \ge 0 \end{array}\right.$.

Thus, we have one-side Fourier Transformation

$$ \hat F(\omega) = \int_{-\infty}^{\infty} F(t) e^{-i\omega t}dt = \int_0^{\infty}f(t)e^{-\gamma t}e^{-i\omega t}dt = \int_0^{\infty}f(t)e^{-(\gamma +i\omega) t}dt = \int_0^{\infty}f(t)e^{-st}dt = \bar f(s) \tag{32}$$

On the contrast,

$$F(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \hat F(\omega)e^{i\omega t}d\omega $$

So we can get,

$$f(t) = e^{\gamma t}F(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \bar f(s)e^{(\gamma+i\omega) t}d\omega = \frac{1}{2\pi} \int_{-\infty}^{\infty} \bar f(s)e^{st}d\omega =  \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma+i\infty} \bar f(s)e^{s t}ds$$

Means $s = \gamma +i\omega$, $ds = id\omega \to d\omega = \frac{1}{i}ds$.

- We get the Laplace Transform pair,

$$\bar f(s) = \int_0^{\infty}f(t)e^{-st}dt \tag{33} $$
$$f(t) = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma+i\infty} \bar f(s)e^{s t}ds \tag{34}$$

Eqn $(33)$ leads to

$$\cal{L}(df/dt) = s \bar f(s) - f(0) $$
$$\cal{L}(d^2f/dt^2) = s^2 \bar f (s) - s f(0) - f^\prime (0)$$
$$\cal{L}(e^{at}) = \frac{1}{s-a}$$
$$\cal{L}\left\{ f(t) * g(t)\right\} = \bar f(s) \bar g (s)$$
$$\cal{L}\{1\} = \frac{1}{s}$$

So if we have a single degree of freedom system

$$\ddot x + \frac{c}{m} \dot x +\frac{k}{m} x =0 $$

Takes Laplace transform

$$s^2 \bar x -s x(0) -\dot x(0) + \frac{c}{m}s\bar x - x(0) + \frac{k}{m} \bar x = 0$$

If $x(0) = 2$ and $\dot x(0) = -5$, are the initial conditions.

$$(s^2+\frac{c}{m}s+\frac{k}{m}) \bar x(s) = 2s +5$$

The solution is 
$$\bar x (s) = \frac{2s+5}{s^2+\frac{c}{m}s+\frac{k}{m}} $$

- If $$\ddot x + \frac{c}{m} \dot x +\frac{k}{m} x = u(t) $$

then the solution is
$$\bar x (s) = \frac{2s+5}{s^2+\frac{c}{m}s+\frac{k}{m}} + \frac{\bar u}{s^2+\frac{c}{m}s+\frac{k}{m}}$$

