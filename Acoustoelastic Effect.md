@[TOC](目录)

<script type="text/javascript"
src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>

# 声弹性效应(Acoustoelastic effect)

## 声弹性效应(acoustoelastic effect)

是指在初始静态应力场作用下，弹性材料的声速（纵波和横波速度）的变化。这是连续质量材料中机械应力和有限应变之间的非线性本构关系。在经典线弹性理论中，大部分弹性材料的小变形可以用外加应力和产生的应变之间的线性关系来描述。这种关系通常被称为广义胡克定律。线弹性理论涉及到二阶弹性常数（例如，$\lambda$和$\mu$），在弹性材料中产生恒定的纵向和剪切声速，不受外加应力影响。声弹性效应则包括施加应力和产生应变之间的本构关系（非线性弹性理论[1]）的高阶展开式，它产生的纵向和剪切声速依赖于材料的应力状态。在无应力材料的极限下，再现了线弹性理论的声速。

早在1925年，布里渊(Brillouin)就对声弹性效应进行了研究[2]。他发现声波的传播速度会随所施加的静水压力成比例地下降。然而，他的理论的一个结果是声波会在足够大的压力下停止传播。后来证明，这种荒谬的效应是由弹性参数不受压力影响的错误假设造成的[3]。

1937年，Murnaghan[4]提出了一种数学理论，将线弹性理论扩展到包括弹性各向同性材料的有限变形。该理论包括三个三阶弹性常数$l$、$m$和$n$。1953年，Huges and Kelly[5]在他们的实验工作中使用Murnaghan的理论，建立了一些弹性材料的高阶弹性常数的数值，这些材料包括聚苯乙烯、阿姆科铁和派热克斯玻璃，它们受到静水压力和单轴压缩。

## 超弹性(hyperelastic)材料的非线性弹性理论

声弹性效应是非线性弹性材料有限变形的效应。关于这一点的现代全面说明可以在[1]中找到。这本书处理非线性弹性理论的应用、和能够实现大弹性变形的固体材料的机械性能的分析。对于可压缩的各向同性超弹性材料（如多晶钢[6]）的声弹性理论的特殊情况，在本文中根据Ogden[1]提出的非线性弹性理论进行复制和显示。

注意，本文和[1]的设置都是等温条件，没有涉及热力学内容。

### 本构关系-超弹性材料（应力-应变关系）

超弹性材料是一种特殊的Cauchy弹性材料，在这种材料中，任意一点的应力都是客观的，而且只取决于关于任意参考构形的变形的当前状态（关于变形的更多细节，参见【变形（力学）】和【有限应变】）。然而，应力所做的功可能取决于变形的路径。因此，Cauchy弹性材料具有非保守(non-conservative)结构，其应力不能由标量弹性势函数导出。在Cauchy弹性材料的特殊情况下，应力所做的功与变形路径无关，这被称为Green弹性或超弹性材料。这种材料是保守的，材料中的应力可以由标量弹性势导出，通常称为应变能密度函数。

根据所选择的应力和应变形式，应力和应变之间的本构关系可以有不同的表达形式。选择一阶Piola-Kirchhoff应力张量$\bm{P}$（它是名义(nominal)应力张量$\bm{N}=\bm{P}^T$的转置），可压缩超弹性材料的本构方程可以用拉格朗日Green应变（$\bm{E}$）表示：

$\begin{equation} \bm{P} = \bm{F}\cdot \frac{\partial W}{\partial \bm{E}} ~~ \text{or} ~~  P_{ij} = F_{ik}\frac{\partial W}{\partial E_{kj}}\end{equation}$

其中，$\bm{F}$是变形梯度张量(deformation gradient tensor)，而第二个表达式使用爱因斯坦求和约定进行张量的索引表示。$W$是超弹性材料的应变能密度函数，而且已被定义为单位体积而不是单位质量，因为这避免了将等式右边乘以参考构形的质量密度$\rho_0$。[1]

假定在当前应变$\bm{E}$下，标量应变能密度函数$W(\bm{E})$可以被近似为泰勒级数展开式，它可以被表示（用索引表示）为：

$\begin{equation} W = C_0+C_{ij}E_{ij}+\frac{1}{2!}C_{ijkl}E_{ij}E_{kl} + \frac{1}{3!}C_{ijklmn}E_{ij}E_{kl}E_{mn}+ \cdots\end{equation}$

当材料处于未变形状态（即$W(\bm{E}_{ij}=0)=0)$时，施加应变能函数为零且有最小值的限制。可见，应变能函数中既不存在常数项，也不存在线性项，因此：

$\begin{equation} W = \frac{1}{2!}C_{ijkl}E_{ij}E_{kl} + \frac{1}{3!}C_{ijklmn}E_{ij}E_{kl}E_{mn}+ \cdots\end{equation}$

其中，$C_{ijkl}$是二阶弹性模量的四阶张量，而$C_{ijklmn}$是三阶弹性模量的六阶张量。$E_{ij}=E_{ji}$的对称性以及标量应变能密度函数$W$表明，二阶弹性模量$C_{ijkl}$有如下对称性：

$\begin{equation} C_{ijkl} = C_{jikl} = C_{ijlk}\end{equation}$

这将独立弹性常数从81个减少到36个。此外，幂展开表明二阶模量也具有主对称性

$\begin{equation} C_{ijkl} = C_{klij} \end{equation}$

这使独立弹性常数的数目进一步减少到21个。三阶弹性模量$C_{ijklmn}$也可以采用同样的参数。这些对称性也允许用Voigt表示法来表示弹性模量（即，$C_{ijkl} = C_{IJ}$和$C_{ijklmn} = C_{IJK}$）。

变形梯度张量可以用分量形式表示为

$\begin{equation} F_{ij} = \frac{\partial u_{i}}{\partial X_{j}}+\delta_{ji} \end{equation}$

其中，$u_{i}$是材料点$P$从参考构形中坐标$X_i$到变形构形中坐标$x_i$的位移（参见【有限应变理论】页面中的图2）。将应变能函数的幂展开式引入本构关系式，用【有限应变张量】页面给出的展开式取代拉格朗日应变张量$E_{kl}$（注意，在这里用的是小写字母$u$，与【有限应变】页面中的大写字母相比），得到本构方程：

$\begin{equation} P_{ij} = C_{ijkl} \frac{\partial u_{k}}{\partial X_{l}} + \frac{1}{2} M_{ijklmn}\frac{\partial u_{k}}{\partial X_{l}}\frac{\partial u_{m}}{\partial X_{n}} + \frac{1}{3} M_{ijklmnpq}\frac{\partial u_{k}}{\partial X_{l}}\frac{\partial u_{m}}{\partial X_{n}}\frac{\partial u_{p}}{\partial X_{q}} + \cdots \end{equation}$

其中

$\begin{equation} M_{ijklmn} = C_{ijklmn} + C_{ijln}\delta_{km} + C_{jnkl}\delta_{im} + C_{jlmn}\delta_{ik} \end{equation}$

更高阶项被忽略[7][8]（详细推导参见[9]）。通过忽略用$\frac{\partial u _k}{\partial X_l}$表示的高阶项，该表达式退化为$P_{ij} = C_{ijkl} \frac{\partial u_{k}}{\partial X_{l}}$，这是广义胡克定律的一个版本，其中$P_{ij}$是应力的量度，而$\frac{\partial u_{k}}{\partial X_{l}}$是应变的量度，$C_{ijkl}$是二者之间的线性关系。

## 声 速

假设一个小的动态（声学）变形干扰了一个已经受到静态应力作用的材料，声弹性效应可以看作是一个小变形叠加在一个更大的有限变形上的效应（也称为“小-叠-大”理论）[8]。我们定义一个给定材料点的三种状态。在参考（无应力作用）状态中，点由坐标向量$\bm{X}$定义，而在静态应力作用状态下（即，在施加预应力的影响下），同一个点具有坐标向量$\bm{x}$。最后，假设在较小的动态扰动（声应力场）下的材料点有坐标向量$\bm{x}^\prime$。材料点的总位移（在静态预应力和动态声扰动的影响下）可以用位移向量来描述：

$\begin{equation} \bm{u} = \bm{u}^{(0)} + \bm{u}^{(1)} = \bm{x}^\prime - \bm{X} \end{equation}$

其中，$\bm{u}^{(0)}=\bm{x}-\bm{X}$, $\bm{u}^{(1)} = \bm{x}^\prime - \bm{x}$

分别描述了由施加预应力引起的静态（拉格朗日）初始位移和由声扰动引起的（欧拉）位移。对于额外的欧拉扰动的Cauchy第一运动定律（或线性动量平衡），  $\bm{u}^{(1)}$可以用中间拉格朗日变形$\bm{u}^{(0)}$的形式推导出来，假定“小-叠-大”假设

$\begin{equation} \bm{u}^{(0)} >> \bm{u}^{(1)} \end{equation}$

成立。利用Cauchy第一运动定律的拉格朗日形式，其中恒定体力（即重力）的影响被忽略，得到：

注意，下标/上标“0”在本文中使用，以表示无应力参考状态，而带点变量是变量的时间  导数，而  为关于拉格朗日坐标系  的散度算子。

运动定律的等号右端（时间相关部分）可以表示为：


假设无应力状态和初始变形状态都是静态的，因此


对于等式左侧（空间相关部分），关于  的空间拉格朗日偏导数，通过使用连锁规则，可以用欧拉  展开，并且通过位移向量之间的关系式


来改变变量。其中使用了简易形式  。因此


进一步假定静态初始变形  （预应力状态）是处于平衡态的，意味着  ，而运动定律结合上述本构方程，可简化为静态初始变形  和附加动态扰动  之间的线性关系（即  中的高阶项），如[7]（详细推导参见[9]）


其中


这个表达式被认为是线性波动方程。考虑以下形式的平面波：


其中  是传播方向上的拉格朗日单位矢量（即平行于波数  垂直于波前），  为单位向量，被称为偏振向量（描述质点运动方向），  为相速度，  为二次连续可微函数（例如一个正弦函数）。把这个平面波插入到上面推导的线性波动方程中，得到[10]


其中  被引入，作为声张量，取决于  [10]：


这个表达式称为传播条件，它决定了对于给定传播方向  ，对应于平面波的可能波的速度和偏振。波速可由特征方程确定[10]：


其中  为行列式，  为特征矩阵。

对于超弹性材料，  是对称的（但并不是一般情况），而特征值  因此为实数。为了让波速也为实数，特征值需要为正。[1]如果是这样，则在给定的传播方向  上，存在三个相互正交的实平面波。从声张量的两种表达式中可以清楚地看出，对于所有非零向量  ，有[10]：


和不等式  （也被称为强椭圆率条件(strong ellipticity condition)），而且  保证均匀平面波波速为实数。偏振  对应于纵波，其质点运动平行于波传播方向（也被称为压缩波）。两个偏振  时对应于横波，其质点运动垂直于波传播方向（也被称为剪切波）。[10]