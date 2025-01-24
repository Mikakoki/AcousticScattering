# Redirection and Splitting of Sound Waves by a Periodic Chain of Thin Perforated Cylindrical Shells

## Overview
This repository showcases Mathematica code for recoding this [article](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.7.064034). 

## Contents
- Theoretical Analysis
- Mathematical computations
- Plotting examples
- Errors of the Code


## Theoretical Analysis

### Scattering of a Cylinder Object
The first problem we have to solve is the acoustic wave equation scattering from a cylinder.

$$\left(\nabla^2-\frac{1}{c^2} \frac{\partial^2}{\partial t^2}\right) p=0$$

writing the equation in cylindrical coordinates and pertaining the variables to $p(r,\phi,z,t) = R(r) \Phi (\phi) Z(z) T(t)$, we rewrite the wave equation as 

$$\frac 1R \frac{\partial^2 R}{\partial r^2}+\frac{1}{R r} \frac{\partial R}{\partial r}+\frac{1}{r^2 \Phi} \frac{\partial^2 \Phi}{\partial \phi^2}+\frac 1Z \frac{\partial^2 Z}{\partial z^2}	= \frac{1}{c^2 T} \frac{\partial^2 T}{\partial t^2}$$

in which; 

$$\frac 1Z \frac{\partial^2 Z}{\partial z^2}= -k_z^2$$

$$\frac{1}{c^2 T} \frac{\partial^2 T}{\partial t^2} = -k^2$$

$$\frac 1R \frac{\partial^2 R}{\partial r^2}+\frac{1}{R r} \frac{\partial R}{\partial r}+\frac{1}{r^2 \Phi} \frac{\partial^2 \Phi}{\partial \phi^2} = -k^2+k_z^2 = -k_r^2 $$

where $\frac{1}{\Phi} \frac{\partial^2 \Phi}{\partial \phi^2}=- n^2$, therefore for the radial differentiation equation we will have:

$$\frac{\partial^2 R}{\partial r^2}+\frac{1}{ r} \frac{\partial R}{\partial r}+(k_r^2 -\frac{n^2}{r^2})R =0 $$

the solutions for each set of equations are:

$$R(r)=R_1 J_n(k_r r) + R_2 Y_n(k_r r) $$

$$T(t) = T_1 e^{-i\omega t}$$

$$Z(z) = Z_1 e^{ik_z z}$$

$$\Phi(\phi) = \Phi_1 e^{in \phi}$$
where $\omega = kc$, and $n$ is an integer. Multiplying all the equations, for all possible nodes and all possible values for wave number in z axis, to have a solution in two dimensions, we will reach the most general solution:

$$p(r,\phi,z,t) = \sum_{n=-\infty}^{\infty}e^{in\phi}\frac 1{2\pi} \int_{-\infty}^{\infty} dk_z [A_nJ_n(k_r r) + B_nY_n(k_r r)]e^{i k_z z -i \omega t}$$

which is for the standing wave. The other form which is commonly used as well is:

$$ p(r,\phi,z,t) = \sum_{n=-\infty}^{\infty}e^{in\phi}\frac 1{2\pi} \int_{-\infty}^{\infty} dk_z [C_nH^{(1)}_n(k_r r) + D_nH^{(2)}_n(k_r r)]e^{i k_z z -i \omega t}$$

by moving to the frequency domain and the boundary condition for a constant z, we reach that the incident and scattered wave are given as in the next section.

### Scattering From a Chain of Cylinders

```latex
The incoming pressure plane wave \( p(\mathbf{r},t) \) can be found by sovling the wave equation in a cylindrical coordinate. the solution for a plane wave will be found as:

\begin{equation}
\begin{aligned}
p(\mathbf{r}, t) & = p_0 \exp (i \mathbf{k} \cdot \mathbf{r} - i \omega t) = p\left(r_l, \varphi_l, t\right) \\
& = p_0 e^{i k r_l \cos \left(\varphi_l - \theta\right) - i \omega t} = p_0 \sum_{n=-\infty}^{\infty} i^n J_n\left(k r_l\right) e^{i n\left(\varphi_l - \theta\right) - i \omega t}
\end{aligned}
\end{equation}

The factor \( e^{i\omega t} \) is omitted in further calculations. The field scattered by a system of shells is written as a superposition of outgoing cylindrical waves radiated by each scatterer. We write this field in polar coordinates centered at \( x,y = 0 \):

\begin{equation}
p_{\mathrm{sc}}(r, \varphi) = \sum_{l^{\prime}} \sum_{n=-\infty}^{+\infty} B_{l^{\prime} n} H_n\left(k r_{l^{\prime}}\right) e^{i n \varphi_{l^{\prime}}}
\label{2}
\end{equation}

Here, the index \( l' \) numerates the shells in the chain, and \( H_n \) denotes the Hankel function of the first kind. Using Graf's theorem for \( r_l < r_{ll'} \):

\begin{equation*}
H_n\left(k r_{l^{\prime}}\right) e^{in\varphi_{l^{\prime}}} = \sum_{n^{\prime}=-\infty}^{\infty} H_{n-n^{\prime}}\left(k r_{ll'}\right) e^{i (n-n^{\prime})\varphi_{ll'}} J_{n^{\prime}}\left(k r_l\right) e^{i n^{\prime}\varphi_l}
\end{equation*}

Knowing it's a chain of cylinders, we can replace \( r_{ll'} = |l-l'|d \) and \( \varphi_{ll'} = \frac{\pi}{2} \operatorname{sgn}\left(l^{\prime}-l\right) \). Also, by changing \( n' \to -n' \):

\begin{equation*}
H_n\left(k r_{l^{\prime}}\right) e^{in\varphi_{l^{\prime}}} = \sum_{n^{\prime}=-\infty}^{\infty} i^{\left(n+n^{\prime}\right) \operatorname{sgn}\left(l'-l\right)} H_{n+n^{\prime}}\left(k r_{ll'}\right) J_{n^{\prime}}\left(k r_l\right) e^{-i n^{\prime} \varphi_l}
\end{equation*}

Finally, by small adjustments:

\begin{equation}
H_n\left(k r_{l^{\prime}}\right) e^{i n \varphi_{l^{\prime}}} = \sum_{n^{\prime}=-\infty}^{\infty} i^{\left(n+n^{\prime}\right) \operatorname{sgn}\left(l-l^{\prime}\right)} H_{n+n^{\prime}}\left(k\left|l-l^{\prime}\right| d\right) J_{n^{\prime}}\left(k r_l\right) e^{i n^{\prime}\left(\pi - \varphi_l\right)}
\label{3}
\end{equation}

Replacing equation \ref{3} in equation \ref{2}, we will have:

\begin{equation}
\begin{aligned}
p_{\mathrm{sc}}\left(r_l, \varphi_l\right) = & \sum_{n=-\infty}^{+\infty} \Bigg( B_{l n} H_n\left(k r_l\right) e^{i n \varphi_l} + \sum_{l^{\prime} \neq l} B_{l^{\prime} n} \\
& \times \sum_{n^{\prime}=-\infty}^{+\infty} i^{\left(n+n^{\prime}\right) \operatorname{sgn}\left(l-l^{\prime}\right)} H_{n+n^{\prime}}\left(k\left|l-l^{\prime}\right| d\right) J_{n^{\prime}}\left(k r_l\right) e^{i n^{\prime}\left(\pi-\varphi_l\right)} \Bigg)
\end{aligned}
\label{4}
\end{equation}

Pressure inside the \( l \)th cylinder is expanded over the Bessel functions:

\begin{equation}
p_{\mathrm{in}}(r, \varphi) = \sum_{n=-\infty}^{+\infty} C_{ln} J_n\left(k r_{l}\right) e^{i n \varphi_{l}}
\label{5}
\end{equation}

\subsubsection{Boundary Condition}
We know that the impedance of a perforated shell is almost as much as a plane and is given by:

\begin{equation}
Z_p = -\frac{i \omega \rho_0}{\sigma}\left[h + \frac{16 r}{3 \pi}\left(1 - 2.5 \sqrt{\frac{\sigma}{\pi}}\right)\right] 
\label{6}
\end{equation}

Also, from references, we know that:

\begin{equation*}
Z_p = \frac{\Delta p}{v}
\end{equation*}

We are assuming the speed throughout each tube (perforated hole) is constant, so the speed at the beginning and the end is the same. Since impedance is given as the equation above, the boundary conditions are:

\begin{equation}
\left.v_r\right|_{r=a} = \left.v_r\right|_{r=b} = \frac{\left.p\right|_{r=b} - \left.p\right|_{r=a}}{Z_p}
\end{equation}

Here, the radial velocity in the fluid is:

\begin{equation*}
v_r = -\frac{i}{\omega \rho_0} \frac{\partial p}{\partial r}
\end{equation*}

So for the first boundary condition, we will have:

\begin{equation*}
\frac{\partial (p + p_{sc})}{\partial r}\bigg|_{r=a} = \frac{\partial p_{in}}{\partial r}\bigg|_{r=b}
\end{equation*}

Consider: \( A' = \frac{\partial A}{\partial r} \)

\begin{equation*}
\begin{aligned}
p_0 \sum_{n=-\infty}^{\infty} i^n k J'_n\left(k a\right) e^{i n\left(\varphi_l - \theta\right)} &+  
\sum_{n=-\infty}^{+\infty}\left[B_{l n} k H'_n\left(k a\right) e^{i n \varphi_l} + \sum_{l^{\prime} \neq l} B_{l^{\prime} n} \right. \\
&\times \sum_{n^{\prime}=-\infty}^{+\infty} i^{\left(n+n^{\prime}\right) \operatorname{sgn}\left(l-l^{\prime}\right)} H_{n+n^{\prime}}\left(k\left|l-l^{\prime}\right| d\right) k J'_{n^{\prime}}\left(ka\right) e^{i n^{\prime}\left(\pi - \varphi_l\right)} \left. \right] \\
&= \sum_{n=-\infty}^{+\infty} C_{ln} k J'_n\left(k b\right) e^{i n \varphi_{l}}
\end{aligned}
\end{equation*}
to simplify the equation first by switching indices of \( n \) and \( n' \) in the second phrase and then change $n \to -n$ so the exponential go away. \\ 
Notice: $$J_n(x) = (-1)^n J_{-n}(x) \, \, , \, \, \, \, H_n(x) = (-1)^n H_{-n}(x) \, \, , \, \, \, \, \, i^{-in\pi} = (-1)^n$$ 

$$(-1)^{(n-n')}i^{-(n-n')sgn(l-l')} = (-i^{-sgn(l-l')})^{(n-n')} = i^{(n-n')sgn(l-l')} $$

therefore:
\begin{equation}
p_0 i^n e^{-i n \theta} + B_{ln} \frac{H'_n(ka)}{J'_n(ka)} + \sum_{l' \neq l} \sum_{n'=-\infty}^{\infty} B_{l'n'} i^{(n-n')\operatorname{sgn}(l-l^{\prime})} H_{n-n'}(k|l-l'|d) e^{i n \varphi_l} = C_{ln} \frac{J'_n(kb)}{J'_n(ka)}
\label{8}
\end{equation}

The second boundary condition we will have:

\begin{equation*}
\frac{\left.p_{in}\right|_{r=b} - \left.(p + p_{sc})\right|_{r=a}}{Z_p} = -\frac{-i}{\omega \rho_0} \frac{\partial p_{in}}{\partial r}\bigg|_{r=b} 
\end{equation*}

\begin{equation*}
\begin{aligned}
\sum_{n=-\infty}^{+\infty} C_{ln} J_n\left(k b\right) e^{i n \varphi_{l}} - p_0 &\sum_{n=-\infty}^{\infty} i^n J_n\left(k a\right) e^{i n\left(\varphi_l - \theta\right)} -  
\sum_{n=-\infty}^{+\infty}\left[B_{l n} H_n\left(k a\right) e^{i n \varphi_l} + \sum_{l^{\prime} \neq l} B_{l^{\prime} n} \right. \\
&\times \sum_{n^{\prime}=-\infty}^{+\infty} i^{\left(n+n^{\prime}\right) \operatorname{sgn}\left(l-l^{\prime}\right)} H_{n+n^{\prime}}\left(k\left|l-l^{\prime}\right| d\right) J_{n^{\prime}}\left(ka\right) e^{i n^{\prime}\left(\pi - \varphi_l\right)} \left. \right] \\
&= -\frac{ikZ_p}{\omega \rho_0}\sum_{n=-\infty}^{+\infty} C_{ln} J'_n\left(k b\right) e^{i n \varphi_{l}}
\end{aligned}
\end{equation*}

Again by doing the same procedure in the last part we simplify the equation:

\begin{equation}
\begin{aligned}
p_0 i^n e^{-i n \theta} &+ B_{ln} \frac{H_n(ka)}{J_n(ka)} + \sum_{l' \neq l} \sum_{n'=-\infty}^{\infty} B_{l'n'} i^{(n+n')\operatorname{sgn}(l-l')} H_{n+n'}(k|l-l'|d) e^{i n \varphi_l} \\
&= C_{ln} \frac{J_n(kb)}{J_n(ka)} + \frac{iZ_p}{c_0 \rho_0} C_{ln} \frac{J'_n(kb)}{J_n(ka)}
\end{aligned}
\label{9}
\end{equation}
dividing equations \ref{9} and \ref{8}:
\begin{equation*}
 B_{ln} \left(\frac{H'_n(ka)}{J'_n(ka)}-\frac{H_n(ka)}{J_n(ka)}\right) = C_{ln}\left(\frac{J'_n(kb)}{J'_n(ka)}-\frac{J_n(kb)}{J_n(ka)}- \frac{iZ_p}{c_0 \rho_0} \frac{J'_n(ka)}{J_n(ka)}\right)
\end{equation*}
by substituting the equation above in equation \ref{8}, we eliminate $C_{ln}$ from the equations.

\begin{equation*}
\begin{aligned}
p_0 i^n e^{-i n \theta} &+ B_{ln} \frac{H'_n(ka)}{J'_n(ka)} + \sum_{l' \neq l} \sum_{n'=-\infty}^{\infty} B_{l'n'} i^{(n+n')\operatorname{sgn}(l-l^{\prime})} H_{n-n'}(k|l-l'|d) e^{i n \varphi_l} \notag \\
&= - B_{ln} \frac{\frac{H'_n(ka)}{J'_n(ka)} - \frac{H_n(ka)}{J_n(ka)}}{\frac{J'_n(kb)}{J'_n(ka)} - \frac{J_n(kb)}{J_n(ka)} - \frac{iZ_p}{c_0 \rho_0} \frac{J'_n(ka)}{J_n(ka)}}
\end{aligned}
\end{equation*}
which can be rewritten as:
\begin{equation}
\mathcal{S}_n B_{l n}+\sum_{l^{\prime} \neq l} \sum_{n^{\prime}=-\infty}^{\infty} i^{\left(n-n^{\prime}\right) \operatorname{sgn}\left(l-l^{\prime}\right)} H_{n-n^{\prime}}\left(k\left|l-l^{\prime}\right| d\right) B_{l^{\prime} n^{\prime}} \quad=-p_0 i^n e^{-i n \theta}
\label{10}
\end{equation}
where:
\begin{equation}
\mathcal{S}_n=\frac{H_n(k a)-H_n^{\prime}(k a)\left(\frac{i Z_p}{\rho_0 c_0}+\frac{J_n(k b)}{J_n^{\prime}(k b)}\right)}{J_n(k a)-J_n^{\prime}(k a)\left(\frac{i Z_p}{\rho_0 c_0}+\frac{J_n(k b)}{J_n^{\prime}(k b)}\right)}
\label{11}
\end{equation}
```


## Example Codes

### Sine Wave Visualization
```mathematica
Plot[Sin[x], {x, 0, 2 Pi}]
```

### Integration Example
```mathematica
Integrate[x^2, x]
```


## References
\bibliography{reference.bib}

