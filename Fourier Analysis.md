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

## 2. Fourier Transform

