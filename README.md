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
results are then compared. Check the file *Assignment_2.py*.

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
The primary objectives were:
To qualitatively capture the evolution of non-dimensional acoustic velocity with time and its dependence
on parameters like Heater power, K, time lag, τ etc. How damping affects the fluctuations is also studied.

• To understand how the energy of the system respond to increasing K and how projecting the velocity onto
multiple modes affect the energy distribution of the system.

• To qualitatively capture the evolution of velocity projected onto various Galerkin modes and how they
respond, i.e, grow or decay by varying certain parameters.

• To obtain the hysteresis through a bifurcation plot of the RMS velocity and Heater power value by
increasing and subsequently decreasing the K value and to observe the Hopf and Fold Bifurcations in the
plot.

Check the file *Assignment_3.py*.

### References:
1) *Kosuhik Balasubramanian and R. I. Sujith*, **Thermoacoustic instability in a Rijke tube: Non-normality
and nonlinearity**, *Physics of Fluids 2.4 (2008)*.
2) *Priya Subramanian and Pankaj Wahi Sathesh Mariappan R. I. Sujith*, **Bifurcation analysis of thermoacoustic instability in a horizontal Rijke tube**,  *International journal of spray and combustion dynamics 2.4
(2010)*.
3) *Xiaochuan Yang, Ali Turan, and Shenghui Lei*, **Thermoacoustic Instability in a Rijke Tube with a Distributed Heat Source**, *Journal Of Thermodynamics (2015)*.



