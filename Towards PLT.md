我们记实数域为 $\mathbb{C}$ ，复数域为 $\mathbb{R}$。记一向量为 $\boldsymbol{\vec{r}}\in \mathbb{R}^n$，它的长度为 $r=|\boldsymbol{\vec{r}}|$，对应的单位向量为 $\boldsymbol{\hat r}=\frac{\boldsymbol{\vec{r}}}{r}$。记它的坐标为 $r_0, r_1,...,r_{n-1}$；在 $\boldsymbol{\vec{r}}$ 为三维向量时，记为 $r_x, r_y,r_z$。

取一复数 $z=a+b\text i$，我们定义它的共轭复数（complex conjugate）为 $z^\star=a-b\text{i}$。

取一复数矩阵（或向量） $\boldsymbol A\in \mathbb C^{n\times n}$。我们定义它的共轭转置 $\boldsymbol A^\dagger$为：

$$
\begin{align}
\boldsymbol A^\dagger=({\boldsymbol A}^\star)^\intercal
\end{align}
$$

即取矩阵（或向量）中所有元素的共轭复数后将它转置。

若矩阵 $\boldsymbol A\in \mathbb C^{n\times n}$为正方形，满足 $\boldsymbol A=\boldsymbol A^\dagger$，且对于所有非零向量 $\boldsymbol{\vec{z}}\in \mathbb{C}^n$都满足 $\boldsymbol{\vec{z}}^\dagger \boldsymbol A \boldsymbol{\vec{z}}>0$，则我们称 $\boldsymbol A$为一正定（positive-definite）矩阵。我们记单位矩阵为 $\boldsymbol I$。

取一矩阵 ${\boldsymbol{\Sigma}}\in \mathbb R^{n\times n}$，我们定义各向异性高斯函数为：

$$
\begin{align}
g^{\boldsymbol{\Sigma}}(\boldsymbol{\vec{r}}) ≜ e^{{1 \above 1pt 2}\boldsymbol{\vec{r}}{\Sigma^{-1}}\boldsymbol{\vec{r}}}
\end{align}
$$

记一函数 $f$的傅里叶变换为：

$$
\begin{align}
\mathscr{F}[f]{}(\boldsymbol{\vec\zeta}) ≜ (\frac {1}{2\pi})^\frac{3}{2} \int_{\mathbb{R}^3}\,d^3\boldsymbol{\vec{r}}\space f(\boldsymbol{\vec{r}})\text{e}^{-\text{i}\boldsymbol{\vec{r}}\cdot\boldsymbol{\vec\zeta}}
\end{align}
$$

我们可以用它定义卷积算子 $*$（即卷积定理）：

$$
\begin{align}
f * h \triangleq \mathscr{F}^{-1}(\mathscr{F}[f]\mathscr{F}[h])
\end{align}
$$

我们记一个参考系为 $[\boldsymbol{\mu}]=(\boldsymbol{\hat x, \hat y, \hat z})$ ，后者为一组标准正交基（也就是满足右手法则且定义一个本地坐标系的三个单位向量）。给定一个单位向量，我们就可以为它构建一个本地参考系。有时，我们会需要正交的变基矩阵（change-of-basis matrices），记为 $\boldsymbol Q_{[\boldsymbol \mu]}\in\mathbb{R}^{3\times 3}$ 。它可以将一个世界坐标系中的坐标转化到本地参考系 $[\boldsymbol \mu]$ 。

我们使用四维的Stokes参数向量来描述光的偏振状态，记为 $\boldsymbol {\vec S}^{[\boldsymbol \mu]}$，以及4 $\times$ 4的Mueller矩阵来描述一个物体对光（以Stokes参数向量描述）的作用，记为 $\boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_o]}}$。该Mueller矩阵可以将一个定义在参考系 $[\boldsymbol {\mu_i}]$的Stokes参数向量转化到 $[\boldsymbol {\mu_o}]$（也就是把光从入射方向转到出射方向），并影响它的偏振状态。

我们记一些常见的Stokes向量为：

$$
\begin{align}
\boldsymbol{\vec{S}_ 0}\triangleq\begin{bmatrix}1 \\ 0 \\ 0 \\ 0\end{bmatrix}, \boldsymbol{\vec{S}_ {\text{LHP}}}\triangleq\begin{bmatrix}1 \\ 1 \\ 0 \\ 0\end{bmatrix}, \boldsymbol{\vec{S}_ {\text{LVP}}}\triangleq\begin{bmatrix}1 \\ -1 \\ 0 \\ 0\end{bmatrix}, \boldsymbol{\vec{S}}_ {\text{c}}(\chi,\varsigma)\triangleq\begin{bmatrix}0 \\ 0 \\ \chi \\ \varsigma\end{bmatrix}
\end{align}
$$

其中， $\chi,\varsigma$为对角偏振度和圆偏振度。

最后，我们定义自相关函数 $r$（Autocorrelation Function）和功率谱密度（Power Spectral Density）$p$为

$$
\begin{align}
r(\boldsymbol{\vec d}) &≜ \int_{\mathbb{R}^3}\,d^3\boldsymbol{\vec{r}}\space f(\boldsymbol{\vec{r}+\vec{d}})f(\boldsymbol{\vec{r}})=
({2\pi})^\frac{3}{2} [f(\boldsymbol{-\vec{r}})*f(\boldsymbol{\vec{r}})](\boldsymbol{\vec d})\\
p(\boldsymbol{\vec\zeta})&\triangleq\mathscr{F}[r]{}(\boldsymbol{\vec\zeta})
\end{align}
$$

我们记一时空中的电场为 $\boldsymbol{\vec{E}}(\boldsymbol{\vec{p}}, t)$。它可以被写为数个单色平面波的叠加，即：

$$
\begin{align}
\boldsymbol{\vec{E}}(\boldsymbol{\vec{p}}, t) = (\frac {1}{2\pi})^\frac{3}{2} \int_{\mathbb{R}^3}\,d^3\boldsymbol{\vec{k}}\space \boldsymbol {\vec {a}}_ \perp(\boldsymbol{\vec{k}})\text{e}^{-\text{i}(\boldsymbol{\vec{k}}\cdot\boldsymbol{\vec{p}}-ckt)}
\end{align}
$$

其中 $k=|\boldsymbol{\vec{k}}|$为波数（wavenumber）； $\boldsymbol {\vec {a}}_ \perp(\boldsymbol{\vec{k}})$ 为一描述平面波的振幅和方向的向量方程，它满足 $\boldsymbol {\vec {a}}_ \perp(\boldsymbol{\vec{k}})\cdot \boldsymbol{\vec{k}}=0$ （即垂直于 $\boldsymbol{\vec{k}}$)。

另记平面波场为 $\boldsymbol{\vec{\Psi}}(\boldsymbol{\vec{p}}, t)=\boldsymbol {\vec {a}}_ \perp(\boldsymbol{\vec{k}})\text{e}^{-\text{i}\omega t}$。其中的每个平面波向波向量（wavevector） $\boldsymbol{\vec{k}}$的方向传播；它的振幅为 $\boldsymbol {\vec {a}}_ \perp(\boldsymbol{\vec{k}})$ ，角频率（angular frequency）为 $\omega=ck$，波长（wavelength）为 $\lambda=\frac{2\pi}{k}$ 。不难看出，电场即为平面波场的傅里叶逆变换。

设想一个照亮场景并产生传播到点 $\boldsymbol{\vec{p}}$的光的光源。我们可以将它当作数个点光源的集合。这些点光源各向点 $\boldsymbol{\vec{p}}$放射一个平面波 $\boldsymbol{\vec{\Psi}}(\boldsymbol{\vec{k}}, t)$（其中 $\boldsymbol{\vec{k}}$自光源指向点 $\boldsymbol{\vec{p}}$）。此时，自光源到达点 $\boldsymbol{\vec{p}}$的电场即为所有这些平面波的叠加。

你或许会想直接以蒙特卡洛积分暴力求解电场。但由于积分中高频的复指数函数，这么做会非常难。

现在令前述的光源被一远大于它的虚数圆（imaginary unit sphere）包围。取该圆上的一小块连续表面 $\Delta S$和它对着（subtend）的立体角 $\Omega$。我们定义辐射入该立体角光谱通量（spectral flux）为：

$$
\begin{align}
\Phi\triangleq\oint_{\Delta S}\,d^2\boldsymbol{\vec{k}}|\vec\Psi(k\boldsymbol{\hat{k}})|^2
\end{align}
$$

（上式其实就是把所有通过 $\Delta S$的平面波的能量积分在一起）

式中波数 $k$和由它定义的角频率 $\omega$和波长 $\lambda$保持不变。由于 $|\text{e}^{-\text{i}\omega t}|=1$，我们可以略去原式中的时间。

我们称所有通过 $\Delta S$（或辐射入 $\Omega$）的平面波的集合为一个波袋（wave packet）或光束（beam）（两者可以互相替换）。

波袋拥有清晰的传播方向和不变的波长（以及 $\text{e}^{-\text{i}\omega t}$）。因此，在使用波袋描述光时，我们只需考虑它的波长，而无需考虑时间和时域相干性。追踪固定波长的波袋的传播（就像光线追踪那样），而非一个基于时间的场（像是电场）的演变，这也更符合我们对渲染的认知。



现在，我们将把辐射度量学中的irradiance和radiance泛化（generalize）到两点形式（two-point formalism）以准确地描述波袋/光束的统计数据。

记整个光束的平均传播方向为单位向量 $\boldsymbol{\hat r}$，并以它构建参考系 $(\boldsymbol{\hat x, \hat y, \hat r})$。

我们定义电场的两个横向分量为 $E_x(\boldsymbol{\vec{p}},t)=\boldsymbol{\vec{E}}(\boldsymbol{\vec{p}}, t)\cdot \boldsymbol{\hat x}$和 $E_y(\boldsymbol{\vec{p}},t)=\boldsymbol{\vec{E}}(\boldsymbol{\vec{p}}, t)\cdot \boldsymbol{\hat y}$，并定义它们间的四个互谱密度函数 (cross-spectral density, CSD)为：

$$
\begin{align}
C_{\alpha\beta}(\boldsymbol{\vec{p}},\boldsymbol{\vec{\xi}};\omega)\triangleq\langle{E_\alpha(\boldsymbol{\vec{p}}+\frac{1}{2}\boldsymbol{\vec{\xi}})E_\beta(\boldsymbol{\vec{p}}+\frac{1}{2}\boldsymbol{\vec{\xi}})^\star}\rangle_\omega
\end{align}
$$

其中 $\alpha,\beta \in \{x,y\}$。算子 $\langle\cdot\rangle_\omega$代表系综平均（即在时间上的平均/期望值，式中因此省略了时间），其中下标 $\omega$表示我们在对波包中的所有频率一样的平面波作平均。（日常尺度上电场可以当作一个完全随机的过程，因此我们实际观测到的都是时间上的平均值）

CSD函数可以描述到达两点 $\boldsymbol{\vec{p}}\pm \frac{1}{2}\boldsymbol{\vec{\xi}}$ 且波长为 $\omega$的波的相似度。其中 $\boldsymbol{\vec{p}}$为空间中的任意一点， $\boldsymbol{\vec{\xi}}$则定义两点的位置之差。

有了CSD函数，我们便可以前面提到的Stokes参数向量泛化到两点形式。我们定义泛化Stokes参数（generalized Stokes parameters, gSP）向量为：

$$
\begin{align}
\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu]}}}(\boldsymbol{\vec{p}},\boldsymbol{\vec{\xi}};\omega)
\triangleq\begin{bmatrix}C_{xx} + C_{yy}\\ C_{xx} - C_{yy}\\ C_{yx} + C_{xy}\\ \text{i}(C_{yx} + C_{xy})\end{bmatrix}
\end{align}
$$

上标 $\leftrightarrow$代表泛化后的量。gSP向量与irradiance单位相同且为后者的两点形式泛化（因此它也可以称作泛化irradiance）。它能方便地表示光束的各种数据，比如它的相干性（coherence），强度（intensity），和偏振（polarization）。它的第一项 $S_0^{[\boldsymbol{\mu}]}$为波袋的强度，而它的后三项则描述波袋的偏振状态。同时，该向量的四个项还满足：

$$
\begin{align}
S_0^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{p}},0)^2\geq S_1^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{p}},0)^2 + S_2^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{p}},0)^2 + S_3^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{p}},0)^2
\end{align}
$$

其中 $S_n^{[\boldsymbol{\mu}]}$代表 $\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu]}}}$的第n项。

在 $\boldsymbol{\vec{\xi}}$变化时，gSP会变化的很快。但在 $\boldsymbol{\vec{p}}$变化时，gSP的变化则会很慢。 因此，我们可以假设光束中的波都向它们的平均方向传播，即 $\boldsymbol{\hat{p}}\approx\boldsymbol{\hat{r}}$，并估计gSP为：

$$
\begin{align}
\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu]}}}(\boldsymbol{\vec{p}},\boldsymbol{\vec{\xi}})\approx\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu]}}}(\boldsymbol{\vec{r}},\boldsymbol{\vec{\xi}})
\end{align}
$$

由前式，我们定义泛化radiance为：

$$
\begin{align}
\begin{split}
\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{r}},\boldsymbol{\vec{\xi}};\omega)&\triangleq
\frac{\partial^2}{\partial A\partial\Omega}\oint_{r\Delta S}\,d^2\boldsymbol{\vec{p}}\space\mathcal{\boldsymbol{\overleftrightarrow{S}}}^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{p}},\boldsymbol{\vec{\xi}};\omega)\\
&\approx r^2\frac{\partial}{\partial A}\mathcal{\boldsymbol{\overleftrightarrow{S}}}^{[\boldsymbol{\mu}]}(\boldsymbol{\vec{r}},\boldsymbol{\vec{\xi}};\omega)
\end{split}
\end{align}
$$

和gSP相似， $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$也拥有与radiance相同的单位（在表面上积分irradiance $\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu]}}}$即为 $\Phi$(spectral flux)，因此 $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$即 $\frac{\partial^2\Phi}{\partial A \partial \Omega}$；现在它看着是不是眼熟多了？）。显然，在 $\boldsymbol{\vec{\xi}}=0$时， $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$的第一项 $\mathcal{{L}}_ 0^{[\boldsymbol{\mu}]}$即为经典的（也就是大家都学过的）radiance $L$（gSP同理）。就像gSP将辐射度量学中的irradiance泛化到两点形式， $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$也将radiance泛化到了两点形式。

$\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$还有两个与 $L$很相似的性质：

性质1： $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$在传播过程中守恒。

性质2：取一对泛化radiance $\boldsymbol{\mathcal{\overleftrightarrow{L}}}_ 1^{[\boldsymbol{\mu}]}$和 $\boldsymbol{\mathcal{\overleftrightarrow{L}}}_ 2^{[\boldsymbol{\mu}]}$。若它们互不相干，则它们的叠加是线性（linear）的（也就是可以直接加在一起）。

$\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$虽然可以完整的描述一个光束，但它的每一项都是方程，而非像辐射度量学那样的具体数值。因此，我们使用各向异性高斯函数展开CSD方程，并以它估计光束的相干性，将 $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$写为：

$$
\begin{align}
\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}(r,\boldsymbol{\vec{\xi}};\omega)&\triangleq
\text e^{\text i k\xi_z}[L_x g^{\Theta_x}(\boldsymbol{\vec{\xi}})\boldsymbol{\vec{S}_{\text{LHP}}}+
L_x g^{\Theta_y}(\boldsymbol{\vec{\xi}})\boldsymbol{\vec{S}_{\text{LVP}}}+
\sqrt{L_xL_y} g^{\Theta_\frac{1}{2}}(\boldsymbol{\vec{\xi}})\boldsymbol{\vec{S}}_{\text{c}}(\chi,\varsigma)]
\end{align}
$$

其中 $r$为平均传播距离， $g$为各向异性高斯函数， $\boldsymbol{\vec{S}}$为常见偏振状态（它们的定义见第一章）。每个横向分量承载的spectral radiance为 $L_x$和 $L_y$，而 $\chi,\varsigma$各为对角和圆偏振度。 $[\boldsymbol{\mu}]=(\boldsymbol{\hat x, \hat y, \hat k})$为本地参考系，其中 $\boldsymbol {\hat k}$（即参考系的正z方向）为波袋的平均传播方向。波袋的相干性由形状矩阵（shape matrices）定义：

$$
\begin{align}
\Theta_{x, y, \frac {1}{2}}\triangleq\begin{bmatrix}\frac{r^2}{k^2}\theta_{x, y, \frac {1}{2}} & \space \\ \space & \sigma_{zz}^2\end{bmatrix}
\space\space\space\text{with}\space\space\space
\theta_\frac{1}{2}\triangleq\frac{\theta_x+\theta_y}{2}
\end{align}
$$

其中 $\theta_{x,y}$各为一正定实数 2$\times$2矩阵。

该光束承载的总spectral radiance即为 $L_x+L_y$，而 $L_x-L_y$， $\chi$，和 $\varsigma$则将光束的能量分配到不同到偏振状态上。在波长一定的情况下，这些数值都为常量。

矩阵 $\Theta_{x, y}$可以描述光束的空间相干性。它的二阶顺序主子式（即左上角 $2\times 2$） $\theta_{x,y}$描述相干面积（coherence area），而 $\sigma_{zz}$则描述相干长度（coherence length），见下图：

fig 3

在光束传播时，这个矩阵也保持不变。

使用上述形式代表波袋时，我们只需要记录 $L_x$， $L_y$， $\chi$， $\varsigma$，以及形状矩阵的变化，但仍旧可以很好的描述光的相干性和偏振性。

此外，波袋的相干面积 $|\theta_{x,y}|$即使在与其他物体交互（如反射，散射）时也会守恒（证明见原文Ch 3.4）。

我们取一泛化irradiance $\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}$并令它照射一个散射体（scatterer，即能让光发生散射的物体）。记它的形状矩阵，横向分量，和两个偏振度为 $\Theta_{x, y, \frac {1}{2}}^{(i)},S_{x, y}^{(i)},\chi^{(i)}$和 $\varsigma^{(i)}$。记（平均）入射方向为 $\boldsymbol{\hat s}$，散射（出射）方向为 $\boldsymbol{\hat r}$，并依这两个方向构建入射和散射本地参考系为 $[\boldsymbol{\mu_i}]=(\boldsymbol{\hat x, \hat y, \hat s})$和 $[\boldsymbol{\mu_o}]=(\boldsymbol{\hat x, \hat y, \hat r})$。

若要描述一个物体的散射性质，我们需要它的PSD $p(\boldsymbol{\vec\zeta};\omega)$和偏振BSDF $\boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_o]}}(\boldsymbol{\hat s},\boldsymbol{\hat r};\omega)$，后者为一Mueller矩阵。

我们记：

$$
\begin{align}
\Xi_{x,y,\frac{1}{2}}^{(i)}\triangleq
\boldsymbol{Q}_{[\boldsymbol{\mu_i}]}
\boldsymbol{\Theta}_{x,y,\frac{1}{2}}^{(i)}
\boldsymbol{Q}_{[\boldsymbol{\mu_i}]}^\intercal
\space\space\space\text{and}\space\space\space
\boldsymbol{\vec h}\triangleq k(\boldsymbol{\hat r} + \boldsymbol{\hat s})
\end{align}
$$

其中 $\boldsymbol{Q}_ {[\boldsymbol{\mu_i}]}$矩阵将形状矩阵转换到物体的本地参考系，并定义散射算子 $\mathscr{D}$为：

$$
\begin{align}
\mathscr{D}(\boldsymbol \Sigma)\triangleq\bigg(\frac {|\Sigma|}{8\pi^3}\bigg)^\frac {1}{2}\Big(\mathcal p*g^{\boldsymbol{\Sigma}^{-1}}\Big)\Big(\boldsymbol{\vec h}\Big)
\end{align}
$$

即物体的PSD和一矩阵为 $\boldsymbol \Sigma$的高斯函数的卷积。

我们定义散射波BSDF（wBSDF）为：

$$
\begin{align}
\begin{split}
\boldsymbol{\mathcal{W}_\text{scat}}
\bigg(\mathcal{\boldsymbol{\overleftrightarrow{S}^{[\mu_i]}}}\bigg)&\triangleq
\cos\vartheta_o \boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_o]}} \\&\times 
\Big[S_x^{(i)} \mathscr{D}\{\boldsymbol \Xi_x^{(i)}\}\boldsymbol{\vec{S}_{\text{LHP}}}+
S_x^{(i)} \mathscr{D}\{\boldsymbol \Xi_y^{(i)}\}\boldsymbol{\vec{S}_{\text{LVP}}}+
\sqrt{S_x^{(i)}S_y^{(i)}} \mathscr{D}\{\boldsymbol \Xi_{\frac {1}{2}}^{(i)}\}\boldsymbol{\vec{S}}_{\text{c}}(\chi,\varsigma)\Big]
\end{split}
\end{align}
$$

其中 $\vartheta_o$为出射倾角（inclination angle）。正如经典BSDF将irradiance转化为radiance，wBSDF也可以将泛化irradiance转为散射radiance。wBSDF可以很好的表示相干性导致的光谱和偏振变化。wBSDF输出的结果为一承载radiance和偏振状态的的经典Stokes参数向量 $\boldsymbol{\vec{L}}^{[\boldsymbol{\mu_0}]}$（即 $\boldsymbol{\vec \xi}=0$时的 $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$），它满足：

$$
\begin{align}
L_{x,y}^{(o)}=\frac{L_0^{[\boldsymbol \mu _o]}\pm L_1^{[\boldsymbol \mu _o]}}{2},\space\space
\chi^{(o)}=\frac{L_2^{[\boldsymbol \mu _o]}}{\sqrt{L_{x}^{(o)}L_{y}^{(o)}}},\space\space
\varsigma^{(o)}=\frac{L_3^{[\boldsymbol \mu _o]}}{\sqrt{L_{x}^{(o)}L_{y}^{(o)}}}
\end{align}
$$

（仔细看看 $\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu}]}$的高斯函数定义你就知道这个关系是怎么出来的了）。

下面，我们讨论散射过程对形状矩阵的影响。

根据Van Cittert-Zernike 定理，一个光束在横截面上的强度分布（intensity distribution）即为它的空间相干性的傅里叶变换，反之亦然。换言之，若一光束的强度分布为各向同性，那么它的空间相干性也应为各向同性，即 $\theta_{x,y}$应正比于 $I$。但若散射过程使强度分布产生各向异性，则透过傅里叶变换，空间相干性会产生与之相反的各向异性（即拉伸变为收缩，反之亦然）。

由于相干面积 $|\theta_{x,y}|$守恒，散射过程只能改变它的离心率（椭圆）和它在横截面上的方向。这两种变形都是几何变形。因此，我们可以通过散射强度的曲率（curvature，也就是二阶导数），即pBSDF的二阶导数，来推导强度分布的变形。它的逆则为相干面积的变形。

对于一个Mueller矩阵，我们可以使用Lu-Chipman分解将它拆成：

$$
\begin{align}
\boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_o]}}=m_{00}\boldsymbol M_{\Delta p}\boldsymbol M_d \boldsymbol M_r
\end{align}
$$

其中 $m_{00}$为平均强度系数，也是矩阵左上角的元素。 $\boldsymbol M_{\Delta p}$为归一化的去偏振器（depolarizer）。 $\boldsymbol M_d$为归一化的偏振双向衰减器（diattenuator）。 $\boldsymbol M_r$为相位延迟器（retarder）。

我们记 $\boldsymbol{\hat r}^\perp$为一垂直于 $\boldsymbol{\hat r}$的面（即光束的横截面）。此时，在横截面上的散射强度的变换即为偏振双向衰减器的黑塞矩阵（Hessian，一个由多元函数的二阶偏导数构成的矩阵）：

$$
\begin{align}
\begin{split}
&\boldsymbol U_{x,y}\triangleq \frac{1}{|{\boldsymbol{\tilde{U}}_ {x,y}^{-1}}|}{\boldsymbol{\tilde{U}}_ {x,y}^{-1}}\space,\space \text{with} \\
&{\boldsymbol{\tilde{U}}_ {x,y}}\triangleq 
\boldsymbol {\vec S}_ 0^{\intercal} \frac{\text d}{\text d \boldsymbol{\vec r '}^2}
\bigg[\cos \vartheta_0 m_{00}\boldsymbol{M}_ d \bigg]\Bigg|_ {\boldsymbol{\hat s},\boldsymbol{\hat r}+\boldsymbol{\vec r '}} 
\boldsymbol{\vec S}_ \text{LHP,LVP}
\end{split}
\end{align}
$$

式中对 $\boldsymbol{\vec r '}\in\boldsymbol{\hat r}^\perp$求导且令 $\boldsymbol{\hat s}$不变。最右边的 $\boldsymbol{\vec S}_ \text{LHP,LVP}$分别对应x和y两个横向分量，而 $\boldsymbol {\vec S}_ 0^{\intercal}$则屏蔽掉式中的总散射radiance（因为它和强度/偏振分布无关）。对 ${\boldsymbol{\tilde{U}}_ {x,y}}$取逆对应前述的傅里叶关系，而最后的归一化则保证相干面积不变。

我们记 $\varphi_z$为相位延迟器 $\boldsymbol M_r$描述的关于z轴的旋转角度。它可以描述整个pBSDF和相位延迟所导致的角度变化。现在，我们可以将散射形状矩阵写为：

$$
\begin{align}
\boldsymbol \theta_{x,y}^{(o)}&=\boldsymbol{U}_ {x,y} \boldsymbol{R}(\varphi_z)^\intercal
(\cos^2\varphi_z \space \boldsymbol \theta_{x,y}^{(i)}+\sin^2\varphi_z \space \boldsymbol \theta_{y,x}^{(i)})
\boldsymbol{R}(\varphi_z) \boldsymbol{U}_ {x,y}^\intercal\\
\boldsymbol R_T&\triangleq\begin{bmatrix}
\cos 2\varphi_z & -\sin 2\varphi_z \\ 
\sin 2\varphi_z & \cos 2\varphi_z\end{bmatrix}
\end{align}
$$

$\boldsymbol \theta_{x,y}^{(o)}$仍为一正定矩阵。上述等式首先基于入射和散射横截面间的旋转将两个入射形状矩阵（ $\boldsymbol \theta_{x}^{(i)}$和$\boldsymbol \theta_{y}^{(i)}$)混合，再将它们转换到出射（散射）方向的横截面上，并最后施加前述的几何变形。

相对roughness和specular lobe

现代散射理论将自散射体出射的能量来源划分为两种：一种是前面提到的散射，而另一种则来自于直接，或高光（specular）反射。后者产生的radiance为：

$$
\begin{align}
\boldsymbol{\mathcal{W}}_ \text{direct}
\Big (\boldsymbol {\vec S}^{(i)}\Big )\triangleq
\boldsymbol M_{\text direct}^{\boldsymbol{[\mu_i]\to[\mu_o]}}\boldsymbol {\vec S}^{(i)} 
\end{align}
$$

其中 $\boldsymbol {\vec S}^{(i)}$为经典irradiance（或Stokes参数向量），即 $\boldsymbol{\vec \xi}=0$时的 $\boldsymbol{\mathcal{\overleftrightarrow{S}}}^{[\boldsymbol{\mu}]}$。pBSDF $\boldsymbol M_{\text direct}^{\boldsymbol{[\mu_i]\to[\mu_o]}}$即为将高光反/折射中的菲涅尔关系写为Mueller矩阵。同时，由于高光反射不会导致各向异性， $\boldsymbol{U}_ {x,y}\equiv \boldsymbol I$。

高光反射和散射的比例取决于物体本身的粗糙度（roughness）。我们定义rms roughness $q$和相对roughness $q_{rel}$为：

$$
\begin{align}
q&\triangleq\bigg(\int_{\mathbb{C}^3}\,\text d^3\boldsymbol{\vec r}\space f(\boldsymbol{\vec r})^2\bigg)^{\frac{1}{2}}=\sqrt{r(0)}\\
q_{rel}&\triangleq\bigg(\frac{1}{2\pi}\bigg)^{\frac{3}{4}}\bigg(\int_{\mathcal{S}^2}\,\text d^2\boldsymbol{\hat r}\space p(k\boldsymbol{\hat s}+k\boldsymbol{\hat r})\bigg)^{\frac{1}{2}}
\end{align}
$$

将两种能量来源结合，我们可以定义总wBSDF为：

$$
\begin{align}
\boldsymbol{\mathcal{W}}&\triangleq
\alpha_\text{direct} \boldsymbol{\mathcal{W}}_\text{direct}\Big(\boldsymbol {\vec S}^{(i)}\Big)+
(1-\alpha_\text{direct} )\boldsymbol{\mathcal{W}}_ \text{scat}
\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}\Big) \\
\alpha_\text{direct}&=\text{exp}[-(k \space\cos(\vartheta_i)\space q_\text{rel})^2]\\
\end{align}
$$

其中 $k$为先前提到的波数； $p$为PSD； $\vartheta_i$为入射倾角。

wBSDF的重要性采样

重要性采样是光线追踪中不可或缺的一环。但当我们自镜头向场景中追踪光束时，我们并不知道它在某一点的形状矩阵，也因此无法进行完全的重要性采样。


The scattering wBSDF rewritten. 

$$
\begin{align}
\boldsymbol{\mathcal{W}}_ \text{scat}
\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}\Big)\triangleq
\boldsymbol{\mathcal{W}}_ \text{MInc}
\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}\Big)+
\boldsymbol{\widetilde{\mathcal W}}_ \text{scat}
\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}\Big)
\end{align}
$$

Maximally-incoherent wBSDF

$$
\begin{align}
\boldsymbol{\mathcal{W}}_ \text{MInc}
\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}^{[\mu_i]}}}\Big) \triangleq
\boldsymbol{\mathcal{W}}_ \text{scat}
\Bigg |_{{\Xi_{x,y,\frac{1}{2}}^{(i)}}\to{\sigma_\text{min} \boldsymbol{I}}}
\end{align}
$$

Shape Matrix of Spontaneous emission sources

$$
\begin{align}
\boldsymbol{\Theta}_ {x,y}=
\begin{bmatrix}
\frac{r^2}{k^2}\frac{2\pi\Omega}{A} & \space & \space\\
\space & \frac{r^2}{k^2}\frac{2\pi\Omega}{A} & \space \\
\space & \space & \sigma_{zz}^2
\end{bmatrix}
\end{align}
$$


Modified transport equation, 其中 $\mathfrak{S}_+^2$ 表示上半球

$$
\begin{align}
\begin{split}
\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu_o}]} &=
\int_{\mathfrak{S}_ +^2}\,\text{d}^2\hat{\boldsymbol{s}}\space\boldsymbol{\mathcal{W}}\Big
(\boldsymbol{\mathcal{\overleftrightarrow{L}}}^{[\boldsymbol{\mu_o}]}\Big)\space
\hat{\boldsymbol{s}}\cdot\hat{\boldsymbol{n}}\\
&=\int_{\mathfrak{S}_ +^2}\,\text{d}^2\boldsymbol{\mathcal{W}}\Big(\boldsymbol{\mathcal{\overleftrightarrow{S}}}^{[\boldsymbol{\mu_o}]}\Big)\space
\end{split}
\end{align}
$$


Gaussian PSD

$$
\begin{align}
p(\boldsymbol{\vec\zeta})\triangleq q^2 \sigma_s^2 g^{\sigma_s^{-2}}(\boldsymbol{\vec\zeta})
\end{align}
$$

Diffraction operator 

$$
\begin{align}
\mathcal{D}\{\boldsymbol{\Sigma}\}=\frac{q^2}{2\pi}{|\boldsymbol{\Sigma}^{-1}+\sigma_s^{-2}\boldsymbol I|^{-\frac{1}{2}}}g^{\boldsymbol{\Sigma}^{-1}+\sigma_s^{-2}\boldsymbol I}\Big(\boldsymbol{\vec h}\Big)
\end{align}
$$

Scattering pBSDF Matrix, in which $f_{pp}$ are the Felipe operator 

$$
\begin{align}
\begin{split}
\boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_{sp}]}}&=\frac {1}{\lambda^2}
\begin{bmatrix}
m_{00} & m_{01} & \space &\space \\
m_{01} & m_{00} & \space &\space \\
\space & \space & m_{22} &m_{23} \\
\space & \space & -m_{23} &m_{22} \\
\end{bmatrix}
\boldsymbol T^{\boldsymbol{[\mu_i]\to[\mu_{sp}]}} \\
\text{with}\space\space m_{00}&=\frac{|f_{pp}|^2+|f_{ss}|^2}{2}, \space\space\space \space m_{01}=\frac{|f_{pp}|^2+|f_{ss}|^2}{2}\\
m_{22}&=\text{Re}\{f_{ss}f_{pp}^\star\}\space\space and\space\space\space m_{23}=\text{Re}\{f_{ss}f_{pp}^\star\}
\end{split}
\end{align}
$$

Hessians of the pBSDF

$$
\begin{align}
{\boldsymbol{\tilde{U}}_ {x}}= 
\frac{\text d}{\text d \boldsymbol{\vec r '}^2}
\Big|f_{pp}\Big|^2\Bigg|_ {\boldsymbol{\hat s},\boldsymbol{\hat r}+\boldsymbol{\vec r '}}
\space\space\space \text and \space\space\space
{\boldsymbol{\tilde{U}}_ {y}}= 
\frac{\text d}{\text d \boldsymbol{\vec r '}^2}
\Big|f_{ss}\Big|^2\Bigg|_ {\boldsymbol{\hat s},\boldsymbol{\hat r}+\boldsymbol{\vec r '}}
\end{align}
$$

Importance sample-Diffraction operator

$$
\begin{align}
\mathcal{D}\{\sigma_{\text {min}} \boldsymbol{I}\}=
\frac{q^2}{2\pi w}g^w(h_x)g^w(h_y),\space\space\text{with} \space\space w=\sigma_{min}^{-1}+\sigma_{s}^{-2}
\end{align}
$$

Perfect Diffuse pBSDF Matrix

$$
\begin{align}
\boldsymbol M^{\boldsymbol{[\mu_i]\to[\mu_{sp}]}}=\frac {m_{00}}{\pi\lambda^2}\boldsymbol T^{\boldsymbol{[\mu_i]\to[\mu_{sp}]}}
\end{align}
$$

Complex Transmission Function for Diffraction Grating

$$
\begin{align}
\phi(x)\triangleq \text{exp}\space
\bigg[\text{i}k\frac{b}{2\cos{\vartheta_i}\cos{\vartheta_o}}\sin\bigg(2\pi\frac{x}{\Lambda}\bigg)\bigg]
\end{align}
$$

The PSD of Diffraction Grating

$$
\begin{align}
p(\zeta)=\sum_{n=-\infty}^\infty\, J_n
\bigg(k\frac{b}{2\cos{\vartheta_i}\cos{\vartheta_o}}\bigg)^2
\delta\bigg(\zeta - 2\pi\frac {n}{\Lambda}\bigg)
\end{align}
$$

The Known relation

$$
\begin{align}
\sin\vartheta_o+\sin\vartheta_i=\frac{n}{\Lambda}\lambda
\end{align}
$$













