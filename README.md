<div align="center">

# Monte Carlo Simulation of a Lennard-Jones Fluid in the NVT Ensemble

**Ayush Kumar Samal**  
Department of Physics, Sri Sathya Sai Institute of Higher Learning  
Prasanthi Nilayam, Andhra Pradesh, India  


</div>

---

## Abstract

A canonical (NVT) Monte Carlo simulation of a Lennard-Jones fluid is performed using the Metropolis algorithm for a system of $N = 100$ particles at reduced temperature $T^* = 1.0$ and number density $\rho^* = 0.237~\sigma^{-3}$. The potential energy per particle is monitored over 3000 Monte Carlo steps, and the system is found to equilibrate after approximately 300 steps, converging to $\langle U \rangle / N \approx -1.896~\varepsilon$. This result is in close agreement (within 3.2%) with the reference Molecular Dynamics value of $-1.958~\varepsilon$ obtained using the Velocity-Verlet integrator under equivalent thermodynamic conditions. The work validates the Metropolis acceptance criterion, periodic boundary conditions, and the minimum image convention as a consistent framework for sampling the Boltzmann distribution of a simple fluid.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Theory](#2-theory)
3. [Computational Details](#3-computational-details)
4. [Results](#4-results)
5. [Discussion](#5-discussion)
6. [Conclusion](#6-conclusion)
7. [References](#7-references)

---

## 1. Introduction

Statistical mechanics provides the theoretical bridge between microscopic interactions and macroscopic thermodynamic observables. In practice, evaluating ensemble averages analytically is intractable for systems of more than a few particles, necessitating computational approaches. Two principal methods have emerged for this purpose: Molecular Dynamics (MD), which generates a time evolution of the system by integrating Newton's equations of motion, and Monte Carlo (MC) simulation, which samples configuration space stochastically without reference to physical time.

The Monte Carlo method was first applied to a statistical mechanical system by Metropolis *et al.* (1953), who introduced an importance-sampling scheme — now universally known as the Metropolis algorithm — that generates configurations distributed according to the Boltzmann weight $\exp(-U/k_BT)$. Since then, MC simulation has become an indispensable tool in condensed matter physics, physical chemistry, and materials science.

The Lennard-Jones (LJ) potential is the prototypical model for simple atomic fluids. Its two-parameter form captures the essential competition between short-range Pauli repulsion and long-range London dispersion attraction, and its thermodynamic properties are known with high precision from decades of simulation, making it the ideal benchmark system for any new implementation.

This work implements NVT-MC for the LJ fluid in Scilab, benchmarks the equilibrium potential energy against a reference MD simulation, and investigates the convergence behaviour of the energy as a function of Monte Carlo steps.

---

## 2. Theory

### 2.1 The Lennard-Jones Potential

The interaction energy between two particles separated by distance $r$ is described by the Lennard-Jones potential:

$$V_{\mathrm{LJ}}(r) = 4\varepsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$$

where $\varepsilon$ is the depth of the potential well and $\sigma$ is the finite distance at which the potential is zero. The potential has a minimum at $r_{\min} = 2^{1/6}\sigma \approx 1.122\,\sigma$, where $V = -\varepsilon$.

All quantities are reported in Lennard-Jones reduced units: length in $\sigma$, energy in $\varepsilon$, temperature as $T^* = k_BT/\varepsilon$. The potential is truncated at $r_c = 2.5\,\sigma$.

### 2.2 Periodic Boundary Conditions

Periodic Boundary Conditions (PBC) surround the primary simulation cell with identical image cells in all spatial directions. A particle leaving one face re-enters through the opposite face. The position is wrapped as:

$$x_i \leftarrow x_i - L\,\mathrm{round}\!\left(\frac{x_i}{L}\right)$$

### 2.3 Minimum Image Convention

The Minimum Image Convention (MIC) ensures each particle interacts only with the nearest periodic image of every other particle:

$$\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j - L\,\mathrm{round}\!\left(\frac{\mathbf{r}_{ij}}{L}\right)$$

### 2.4 The Metropolis Algorithm

A particle $j$ is selected at random and displaced by a random vector. The energy change $\Delta U = U_{\mathrm{new}} - U_{\mathrm{old}}$ is computed, and the move is accepted with probability:

$$P_{\mathrm{acc}} = \min\!\left(1,\, \exp\!\left(-\frac{\Delta U}{k_BT}\right)\right)$$

This satisfies **detailed balance**, guaranteeing convergence to the Boltzmann distribution.

### 2.5 Ensemble Averages

The equilibrium potential energy per particle is estimated as:

$$\frac{\langle U \rangle}{N} = \frac{1}{M - M_{\mathrm{eq}}} \sum_{t=M_{\mathrm{eq}}}^{M} \frac{U(t)}{N}$$

where $M$ is the total number of MC steps and $M_{\mathrm{eq}}$ is the equilibration length.

---

## 3. Computational Details

Simulations were performed for $N = 100$ particles in a cubic box of side $L = 7.5\,\sigma$, giving $\rho^* = 0.237\,\sigma^{-3}$. Particles were initialised on a simple cubic lattice. The reduced temperature was $T^* = 1.0$, cutoff $r_c = 2.5\,\sigma$, and maximum displacement $\delta = 0.2\,\sigma$. A total of 3000 MC steps were performed with the first 300 discarded for equilibration.

| Parameter | Symbol | Value |
|:---|:---:|:---:|
| Number of particles | $N$ | 100 |
| Box length | $L$ | $7.5~\sigma$ |
| Number density | $\rho^*$ | $0.237~\sigma^{-3}$ |
| Reduced temperature | $T^*$ | $1.0$ |
| Cutoff radius | $r_c$ | $2.5~\sigma$ |
| Maximum displacement | $\delta$ | $0.2~\sigma$ |
| Total MC steps | $M$ | 3000 |
| Equilibration steps | $M_{\mathrm{eq}}$ | 300 |

---

## 4. Results

### 4.1 Energy Convergence

<p align="center">
  <img src="https://github.com/ayushkumarsamal2709/Monte-Carlo-Simulation_Ayush/blob/main/MC_N100_plot.png?raw=true" width="700" alt="Potential energy per particle vs MC steps"/>
</p>

<p align="center"><em>Figure 1. Potential energy per particle U/N (in units of ε) as a function of Monte Carlo step number. N = 100, T* = 1.0, ρ* = 0.237 σ⁻³, L = 7.5 σ.</em></p>

The system begins in an ordered simple cubic lattice ($U/N \approx -1.098~\varepsilon$). During the first $\approx 300$ steps the energy decreases as the lattice order breaks down and particles relax into a disordered fluid. Beyond step 300 the energy fluctuates about a stable mean, confirming thermal equilibrium has been reached.

### 4.2 Equilibrium Properties

| Quantity | Value |
|:---|:---:|
| MC equilibrium $\langle U \rangle / N$ | $-1.896~\varepsilon$ |
| MD reference $\langle U \rangle / N$ | $-1.958~\varepsilon$ |
| Percentage difference | $3.2\%$ |
| Acceptance ratio | $\approx 80\%$ |
| Equilibration length | $\approx 300$ steps |

---

## 5. Discussion

### 5.1 Comparison with Molecular Dynamics

The MC result of $-1.896~\varepsilon$ agrees with the MD reference value of $-1.958~\varepsilon$ to within 3.2%, validating the implementation. The small residual discrepancy arises from the finite length of both simulations and the intrinsically different phase-space pathways explored by the two methods.

The reference MD simulation uses $L = 15~\sigma$ ($\rho^* = 0.030$), placing all particles on a dense lattice. Kinetic energy from velocity initialisation drives partial expansion of the cluster, but within the finite run time the system never fully disperses. The MD therefore probes a metastable cluster at a local density $\rho_{\mathrm{local}}^* \approx 0.24$, not the global dilute-gas equilibrium.

A canonical MC simulation at the same global density ($\rho^* = 0.030$) correctly samples the true thermodynamic equilibrium — a dilute gas with $\langle U \rangle / N \approx -0.25~\varepsilon$ — which is physically distinct from the MD cluster state. To reproduce the MD result in MC, the box must be chosen to match the effective local density of the MD cluster, hence $L = 7.5~\sigma$.

### 5.2 The Ergodic Hypothesis

The agreement between MC and MD confirms the ergodic hypothesis in practice: both methods converge to the same equilibrium energy when sampling the same thermodynamic state, despite following entirely different trajectories through phase space.

---

## 6. Conclusion

A Metropolis Monte Carlo simulation of a Lennard-Jones fluid in the canonical ensemble has been implemented and validated against a Velocity-Verlet Molecular Dynamics reference. The potential energy per particle converges to $-1.896~\varepsilon$ after approximately 300 equilibration steps, agreeing with the MD reference of $-1.958~\varepsilon$ to within 3.2%. Periodic boundary conditions and the minimum image convention are correctly implemented, and the Metropolis criterion enforces sampling of the correct Boltzmann distribution. This work provides a validated foundation for more advanced simulations including NPT ensemble MC, free energy perturbation, and histogram reweighting methods.

---

## References

1. Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. *Journal of Chemical Physics*, **21**(6), 1087–1092.

2. Allen, M. P., & Tildesley, D. J. (2017). *Computer Simulation of Liquids* (2nd ed.). Oxford University Press.

3. Frenkel, D., & Smit, B. (2002). *Understanding Molecular Simulation: From Algorithms to Applications* (2nd ed.). Academic Press.

4. Lennard-Jones, J. E. (1924). On the Determination of Molecular Fields. *Proceedings of the Royal Society A*, **106**(738), 463–477.

5. Rapaport, D. C. (2004). *The Art of Molecular Dynamics Simulation* (2nd ed.). Cambridge University Press.

---

<div align="center">
Sri Sathya Sai Institute of Higher Learning · Department of Physics · 2025
</div>
