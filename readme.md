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
$$
\left(\nabla^2-\frac{1}{c^2} \frac{\partial^2}{\partial t^2}\right) p=0
$$
writing the equation in cylindrical coordinates and pertaining the variables to $p(r,\phi,z,t) = R(r) \Phi (\phi) Z(z) T(t)$, we rewrite the wave equation as 
$$
\frac 1R \frac{\partial^2 R}{\partial r^2}+\frac{1}{R r} \frac{\partial R}{\partial r}+\frac{1}{r^2 \Phi} \frac{\partial^2 \Phi}{\partial \phi^2}+\frac 1Z \frac{\partial^2 Z}{\partial z^2}	= \frac{1}{c^2 T} \frac{\partial^2 T}{\partial t^2}
$$
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
\begin{equation}
p(r,\phi,z,t) = \sum_{n=-\infty}^{\infty}e^{in\phi}\frac 1{2\pi} \int_{-\infty}^{\infty} dk_z [A_nJ_n(k_r r) + B_nY_n(k_r r)]e^{i k_z z -i \omega t}
\end{equation}
which is for the standing wave. The other form which is commonly used as well is:
\begin{equation}
p(r,\phi,z,t) = \sum_{n=-\infty}^{\infty}e^{in\phi}\frac 1{2\pi} \int_{-\infty}^{\infty} dk_z [C_nH^{(1)}_n(k_r r) + D_nH^{(2)}_n(k_r r)]e^{i k_z z -i \omega t}
\end{equation}
by moving to the frequency domain and the boundary condition for a constant z, we reach that the incident and scattered wave are given as in the next section.


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

