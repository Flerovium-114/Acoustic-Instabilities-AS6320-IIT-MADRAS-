# Acoustic-Instabilities-AS6320-IIT-MADRAS-

## Introduction
Part of assignments for the course **AS6320 - Acoustic Instabilities in Aerospace Propulsion**, *guide:* ***Prof. R.I. Sujith***

## Assignment 2:
The objective of this assignment is to analyse the behaviour of one-dimensional acoustic quiescent flow
fields in ducts with a mean temperature gradient. The wave equation for such a case is derived and is
solved using the Runge-Kutta 4th order method. By transforming the derived wave equation to the mean
temperature space, the differential equation is solved analytically using solutions of Bessel’s differential
equation. The solution is obtained by considering a linear temperature profile with various initial conditions.
The solution is used to obtain the relationship of sound propagation in a quarter wave tube with specified
boundary conditions. The variation of acoustic pressure and velocity amplitudes, the location of nodes and
anti-nodes of acoustic pressure and acoustic velocity and so on are plotted. The numerical and analytical
results are then compared. 

The wave equation derived from the linearised continuity, momentum and energy equations is:
```math
    \frac{\partial^{2}p^{'}}{\partial x^{2}} + \frac{1}{\overline{T}}\frac{\partial \overline{T}}{\partial x}\frac{\partial p^{'}}{\partial x} - \frac{1}{\gamma R\overline{T}}\frac{\partial^{2}p^{'}}{\partial t^{2}} = 0
```

Now introducing perturbations of the form:
```math
p^{'}(x,t) = \hat{p}(x)e^{i\omega t}
```
where $$\omega$$ is the complex frequency of the perturbations.
The wave equation can be written as:

```math
\frac{\partial^{2}\hat{p}}{\partial x^{2}} + \frac{1}{\overline{T}}\frac{\partial \overline{T}}{\partial x}\frac{\partial 
 \hat{p}}{\partial x} + \frac{\omega ^{2}\hat{p}}{\gamma R\overline{T}} = 0
```
We can decouple the second order wave equation into two first order ODEs as follows:
```math

    \frac{d\hat{p}}{dx} = z 
    \frac{dz}{dx} = - \frac{1}{\overline{T}}\frac{d\overline{T}}{dx}z - \frac{\omega^{2}}{\gamma R\overline{T}}\hat{p} 
```
Consider a linear variation of temperature:
```math
 T = T_{0} + mx
```
Since the duct is open at one end and closed at the other end, the frequency of oscillations will follow:
```math
f = \frac{nc}{4L}
```
for integer values of $$n$$.
- $$n$$ determines the mode of oscillation.
- $$c$$ is the speed of sound.
- $$L$$ is the length of the duct.
-  $$\omega = 2\pi f$$.
- $$c = \sqrt{\gamma RT}$$.
Now we can solve these using RK-4 integration method.


### References:
1) **AN EXACT SOLUTION FOR ONE-DIMENSIONAL ACOUSTIC FIELDS IN DUCTS WITH AN AXIAL TEMPERATURE GRADIENT**, *R.I. SUJITH, G.A. WALDHERR AND B.T. ZINN*, *Department of Aerospace Engineering, Georgia Institute of Technology, Atlanta, Georgia 30332, USA*.

## Assignment 3:
### Objectives
To qualitatively capture the evolution of non-dimensional acoustic velocity with time and its dependence
on parameters like Heater power, K, time lag, τ etc. How damping affects the fluctuations is also studied.

• To understand how the energy of the system respond to increasing K and how projecting the velocity onto
multiple modes affect the energy distribution of the system.

• To qualitatively capture the evolution of velocity projected onto various Galerkin modes and how they
respond, i.e, grow or decay by varying certain parameters.

• To obtain the hysteresis through a bifurcation plot of the RMS velocity and Heater power value by
increasing and subsequently decreasing the K value and to observe the Hopf and Fold Bifurcations in the
plot.

### Governing equations
The same governing equations are used with a change in the energy equation.
A modified form of King's law is used to model the heat release rate. The empirical model suggested by Heckl is used where using a factor of $$\frac{1}{3}$$ for the $$u_{0}$$ seemed to work.
```math
\dot{\tilde{Q}}' = \frac{2L_{w}(T_{w} - \overline{T})}{S\sqrt{3}}\sqrt{\pi \lambda C_{v}\overline{\rho}\frac{d_{w}}{2}}\left[\sqrt{\left |\frac{u_{0}}{3} + \tilde{u}_{f}'(t - \tau)\right |} - \sqrt{\frac{u_{0}}{3}}\right]\delta_{D}(\tilde{x} - \tilde{x}_{f})
```
where, $$\delta_{D}$$ is the dirac delta function and we can use the following property,
```math
    \delta_{D}(\tilde{x} - \tilde{x}_{f}) = \delta_{D}(L_{a}(x - x_{f})) = \frac{1}{L_{a}}\delta_{D}(x - x_{f})
```
along with $$\tilde{u_{f}} = u_{0}u_{f}'$$.
After adding a damping term, $$\zeta p'$$ and using equation we can rewrite equation as,  
```math
    \frac{\partial{p'}}{\partial{t}} + \gamma M \frac{\partial{u'}}{\partial{x}} + \zeta p' = (\gamma - 1)\frac{2L_{w}(T_{w} - \overline{T})}{S\sqrt{3}c_{0}\overline{p}}\sqrt{\pi \lambda C_{v}\overline{\rho}\frac{d_{w}}{2}u_{0}}\left[\sqrt{\left |\frac{1}{3} + u_{f}'(t - \tau)\right |} - \sqrt{\frac{1}{3}}\right]\delta_{D}(x - x_{f})
```
We can define a constant,
```math
    K' = (\gamma - 1)\frac{2L_{w}(T_{w} - \overline{T})}{S\sqrt{3}c_{0}\overline{p}}\sqrt{\pi \lambda C_{v}\overline{\rho}\frac{d_{w}}{2}u_{0}}
```
and rewrite equation as, 
```math
    \frac{\partial{p'}}{\partial{t}} + \gamma M \frac{\partial{u'}}{\partial{x}} + \zeta p' = K'\left[\sqrt{\left |\frac{1}{3} + u_{f}'(t - \tau)\right |} - \sqrt{\frac{1}{3}}\right]\delta_{D}(x - x_{f}) 
```
$L_{w}$ is the equivalent length of the wire, $$\lambda$$  is the heat conductivity of air, $$C_{v}$$ is the specific heat of air at constant volume, $$\tau$$ is the time lag, $$\rho$$ is the mean density of air, $$d_{w}$$ is the diameter of the wire, $$T_{w} - \overline{T}$$ is the
temperature difference, and $$S$$ is the cross-sectional area of
the duct. $$L_{w} = 2.2 m$$, $$d_{w} = 0.0005m$$, $$S = 1.8 \times 10^{-3} m^{2}$, $L_{a} = 1 m$$, $$u_{0} = 0.5 m/s$$, $$\lambda = 0.0328 W/(m.K)$$, $$\rho = 1.205 Kg/m^{3}$$ as is the case for a typical laboratory Rijke tube.


### References:
1) *Kosuhik Balasubramanian and R. I. Sujith*, **Thermoacoustic instability in a Rijke tube: Non-normality
and nonlinearity**, *Physics of Fluids 2.4 (2008)*.
2) *Priya Subramanian and Pankaj Wahi Sathesh Mariappan R. I. Sujith*, **Bifurcation analysis of thermoacoustic instability in a horizontal Rijke tube**,  *International journal of spray and combustion dynamics 2.4
(2010)*.
3) *Xiaochuan Yang, Ali Turan, and Shenghui Lei*, **Thermoacoustic Instability in a Rijke Tube with a Distributed Heat Source**, *Journal Of Thermodynamics (2015)*.



