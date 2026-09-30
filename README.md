# Two-Class Traffic Flow Model with Creeping

Numerical study of a heterogeneous two-class traffic-flow model based on systems of nonlinear conservation laws.

The project investigates the interaction between **small and large vehicles** in congested traffic, with particular attention to the phenomenon of **creeping**: small vehicles, such as motorcycles, may continue moving through available spaces even when larger vehicles are completely stopped.

The model is based on the heterogeneous multiclass traffic model proposed by **Fan and Work (2015)** and is compared with the classical **n-population traffic-flow model**.

## Project Overview

Classical macroscopic traffic models, such as the **Lighthill–Whitham–Richards (LWR) model**, describe traffic as a continuous flow governed by a conservation law:

\[
\partial_t \rho + \partial_x(\rho V(\rho)) = 0,
\]

where:

- \(\rho(x,t)\) is the traffic density;
- \(V(\rho)\) is the velocity-density relation;
- \(Q(\rho)=\rho V(\rho)\) is the traffic flow.

The project extends this framework to a heterogeneous traffic system composed of two vehicle classes:

- \(\rho_1\): small vehicles;
- \(\rho_2\): large vehicles.

Unlike standard multiclass models, the two classes are assigned different maximum occupied-space thresholds. This allows the model to reproduce situations in which large vehicles become stationary while smaller vehicles continue moving through the remaining available space.

## Two-Class Creeping Model

The total occupied space is defined as

\[
r = \rho_1 + \rho_2.
\]

The model contains two traffic regimes:

### Non-creeping phase

Both vehicle classes move and satisfy a system of conservation laws:

\[
\begin{cases}
(\rho_1)_t + (\rho_1 V_1(r))_x = 0, \\
(\rho_2)_t + (\rho_2 V_2(r))_x = 0.
\end{cases}
\]

### Creeping phase

When congestion reaches the maximum density allowed for large vehicles, the second class becomes stationary:

\[
\begin{cases}
(\rho_1)_t + (\rho_1 V_1(r))_x = 0, \\
(\rho_2)_t = 0.
\end{cases}
\]

Small vehicles can therefore continue moving even when the large-vehicle population is blocked.

For the numerical experiments, linear Greenshields-type velocity functions are used:

\[
V_1(r)=v_1^{max}\left(1-\frac{r}{r_1^{max}}\right),
\]

\[
V_2(r)=v_2^{max}\left(1-\frac{r}{r_2^{max}}\right).
\]

## Numerical Method

The model is solved numerically using an extension of **Godunov's method** based on the **Cell Transmission Model (CTM)**.

The computational domain is divided into cells and the density of each vehicle class is updated through the numerical fluxes exchanged between neighbouring cells:

\[
\rho_{j,i}^{n+1}
=
\rho_{j,i}^{n}
-
\frac{\Delta t}{\Delta x}
\left(
F_{j,i+1/2}^{n}
-
F_{j,i-1/2}^{n}
\right).
\]

The interface fluxes are determined using **sending and receiving functions**, representing respectively:

- the maximum flow that the upstream cell can send;
- the maximum flow that the downstream cell can receive.

The time step satisfies the CFL stability condition:

\[
v^{max}\frac{\Delta t}{\Delta x} \leq 1.
\]

## Numerical Experiments

Three numerical scenarios are analysed.

### 1. Overtaking

Small vehicles are initially positioned behind large vehicles.

The experiment verifies whether the two vehicle classes can travel at different speeds and whether small vehicles are able to overtake the large-vehicle population.

Both the creeping model and the heterogeneous n-population model reproduce overtaking behaviour.

### 2. Creeping at a Red Traffic Light

Both vehicle populations approach a red traffic light.

Large vehicles accumulate near the boundary and eventually reach their maximum admissible density. In the creeping model, small vehicles can continue moving through the stationary queue.

The comparison with the n-population model highlights the main difference between the two approaches: without a dedicated creeping phase, small vehicles become blocked by the queue of large vehicles.

### 3. Overtaking + Creeping

The final experiment combines the two previous phenomena.

Small vehicles initially overtake part of the large-vehicle population while traffic is moving freely. Once the vehicles approach the red traffic light, large vehicles become stationary and the system transitions from the non-creeping to the creeping phase.

This experiment illustrates how the model can reproduce both **overtaking in free-flow conditions** and **creeping under heavy congestion** within the same mathematical framework.

## Model Parameters

For the creeping model, the numerical experiments use:

\[
v_1^{max}=v_2^{max}=1.8,
\]

\[
r_1^{max}=1.8, \qquad r_2^{max}=1.0.
\]

The computational domain is

\[
x \in [0,50]
\]

with spatial discretisation

\[
\Delta x = 0.05.
\]

The time step is selected according to the CFL condition.

## Main Takeaways

The simulations show that the two-class model can reproduce traffic behaviours that are difficult to represent with standard homogeneous multiclass models:

- interaction between heterogeneous vehicle populations;
- overtaking between vehicle classes;
- different congestion thresholds for small and large vehicles;
- transition between normal traffic and creeping regimes;
- movement of small vehicles through queues of stationary large vehicles.

The comparison with the n-population model demonstrates the importance of using different maximum occupied-space constraints when modelling highly heterogeneous traffic.

## References

- S. Benzoni-Gavage and R. M. Colombo, *An n-populations model for traffic flow*, European Journal of Applied Mathematics, 2003.
- C. F. Daganzo, *The Cell Transmission Model: A Dynamic Representation of Highway Traffic Consistent with the Hydrodynamic Theory*, Transportation Research Part B, 1994.
- S. Fan and D. B. Work, *A Heterogeneous Multiclass Traffic Flow Model with Creeping*, SIAM Journal on Applied Mathematics, 2015.
- S. K. Godunov, *A Difference Scheme for Numerical Computation of Discontinuous Solutions of Hydrodynamic Equations*, 1959.
- M. J. Lighthill and G. B. Whitham, *On Kinematic Waves II: A Theory of Traffic Flow on Long Crowded Roads*, 1955.
- P. I. Richards, *Shock Waves on the Highway*, 1956.
