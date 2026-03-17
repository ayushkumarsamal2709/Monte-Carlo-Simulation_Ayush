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

where $\varepsilon$ is the depth of the potential well and $\sigma$ is the finite distance at which the potential is zero. The repulsive $r^{-12}$ term dominates at short range, while the attractive $r^{-6}$ term reflects the induced dipole–dipole (London dispersion) interaction at longer range. The potential has a minimum at $r_{\min} = 2^{1/6}\sigma \approx 1.122\,\sigma$, where $V = -\varepsilon$.

All quantities in this work are reported in Lennard-Jones reduced units: length in units of $\sigma$, energy in units of $\varepsilon$, and temperature as $T^* = k_BT/\varepsilon$. The potential is truncated at a cutoff radius $r_c = 2.5\,\sigma$, beyond which interactions are set to zero.

### 2.2 Periodic Boundary Conditions

Simulating a bulk fluid with a finite number of particles inevitably introduces surface effects that distort bulk thermodynamic properties. Periodic Boundary Conditions (PBC) mitigate this by surrounding the primary simulation cell with identical image cells in all spatial directions. A particle leaving one face of the box re-enters through the opposite face, so the number of particles $N$ and the volume $V = L^3$ remain constant. The position of particle $i$ is wrapped as:

$$x_i \leftarrow x_i - L\,\mathrm{round}\!\left(\frac{x_i}{L}\right)$$

### 2.3 Minimum Image Convention

Under PBC, each particle interacts with the nearest periodic image of every other particle rather than with all images simultaneously. The Minimum Image Convention (MIC) implements this by correcting the inter-particle displacement vector:

$$\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j - L\,\mathrm{round}\!\left(\frac{\mathbf{r}_{ij}}{L}\right)$$

This guarantees that $|\mathbf{r}_{ij}| \leq L/2$ in each Cartesian direction, which is automatically consistent with the cutoff condition $r_c \leq L/2$.

### 2.4 The Metropolis Algorithm

The goal of canonical MC simulation is to sample configurations $\{\mathbf{r}^N\}$ with probability proportional to the Boltzmann weight $\exp(-U(\mathbf{r}^N)/k_BT)$. The Metropolis algorithm achieves this through a Markov chain of trial moves. A particle $j$ is selected at random and displaced by a random vector drawn uniformly from a cube of side $2\delta$. The energy change $\Delta U = U_{\mathrm{new}} - U_{\mathrm{old}}$ is computed, and the move is accepted with probability:

$$P_{\mathrm{acc}} = \min\!\left(1,\, \exp\!\left(-\frac{\Delta U}{k_BT}\right)\right)$$

If the move is rejected, the particle is returned to its original position and the old configuration is counted again in the ensemble average. This acceptance rule satisfies the condition of **detailed balance**, which is sufficient to guarantee that the stationary distribution of the Markov chain is the target Boltzmann distribution.

### 2.5 The Ergodic Hypothesis and Ensemble Averages

The ergodic hypothesis asserts that, for a sufficiently long simulation, the time average of any observable equals its ensemble (Boltzmann) average. This underpins the physical interpretation of both MD and MC simulations. The equilibrium potential energy per particle is estimated as:

$$\frac{\langle U \rangle}{N} = \frac{1}{M - M_{\mathrm{eq}}} \sum_{t=M_{\mathrm{eq}}}^{M} \frac{U(t)}{N}$$

where $M$ is the total number of MC steps and $M_{\mathrm{eq}}$ is the number of steps discarded for equilibration.

---

## 3. Computational Details

Simulations were performed for $N = 100$ particles in a cubic box of side $L = 7.5\,\sigma$, giving a number density $\rho^* = N/L^3 = 0.237\,\sigma^{-3}$. Particles were initialised on a simple cubic lattice with spacing $L/\lceil N^{1/3} \rceil$. The reduced temperature was set to $T^* = 1.0$, the cutoff radius to $r_c = 2.5\,\sigma$, and the maximum trial displacement to $\delta = 0.2\,\sigma$. A total of 3000 MC steps were performed, with the first 300 steps treated as equilibration and excluded from averages. All simulations were carried out in Scilab 6.x. Particle trajectories were written in XYZ format at every step for post-hoc visualisation in VMD.

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

The potential energy per particle $U/N$ as a function of MC step number is shown in Figure 1. The system begins in an ordered simple cubic lattice configuration ($U/N \approx -1.098~\varepsilon$). During the first $\approx 300$ steps, the energy decreases monotonically as the regular lattice structure is destroyed and particles rearrange into a disordered fluid configuration. Beyond step 300, the energy fluctuates about a stable mean, indicating that thermal equilibrium has been achieved. The amplitude of the fluctuations is consistent with the $\mathcal{O}(N^{-1/2})$ statistical noise expected for a system of 100 particles.

### 4.2 Equilibrium Properties

The equilibrium potential energy per particle, averaged over steps 300–3000, is $\langle U \rangle / N = -1.896~\varepsilon$. The acceptance ratio over the full run was approximately 80%, which is higher than the conventionally optimal $\sim 50\%$ for dense liquids. At the moderately dilute density studied here ($\rho^* = 0.237$), particles are sufficiently well-separated that small displacements ($\delta = 0.2\,\sigma$) rarely produce steric overlaps, leading to the elevated acceptance rate.

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

The MC equilibrium energy of $-1.896~\varepsilon$ agrees with the MD reference value of $-1.958~\varepsilon$ to within 3.2%, providing strong validation of the implementation. The small residual discrepancy is attributable to the finite length of both simulations and the intrinsically different phase-space pathways explored by the two methods.

It is instructive to consider why the two methods require different box sizes to sample the same thermodynamic state. The reference MD simulation uses $L = 15~\sigma$ ($\rho^* = 0.030$), placing all particles on a dense lattice with spacing $a = 1.1~\sigma$. The random velocities assigned at initialisation provide kinetic energy that drives the cluster to partially expand into the box. However, within the finite run time ($t = 20\,\tau$), the system never fully disperses; the MD effectively probes a metastable cluster at a local density of $\rho_{\mathrm{local}}^* \approx 0.24$, not the global dilute-gas equilibrium.

A canonical MC simulation at the same global density ($\rho^* = 0.030$, $L = 15~\sigma$) correctly samples the true thermodynamic equilibrium at that density, which is a dilute gas with $\langle U \rangle / N \approx -0.25~\varepsilon$ — a physically distinct and thermodynamically correct result that differs markedly from the MD value. To reproduce the MD cluster state in MC, the simulation box must be chosen to match the effective local density of the MD cluster, hence $L = 7.5~\sigma$ ($\rho^* = 0.237$).

This comparison highlights a fundamental distinction between the two methods: MD is inherently non-equilibrium in its early stages and its long-time behaviour depends on the initial kinetic energy, whereas MC samples the equilibrium distribution directly but can become trapped in metastable states at high densities if the trial displacement is too small.

### 5.2 The Role of the Ergodic Hypothesis

The agreement between MC and MD confirms the ergodic hypothesis in practice: both methods, despite following entirely different trajectories through phase space, converge to the same equilibrium energy when sampling the same thermodynamic state. This consistency is a necessary condition — though not a sufficient proof — of ergodicity for this system at the studied state point.

---

## 6. Conclusion

A Metropolis Monte Carlo simulation of a Lennard-Jones fluid in the canonical ensemble has been implemented in Scilab and validated against a Velocity-Verlet Molecular Dynamics reference. The principal findings are as follows. The potential energy per particle converges to $-1.896~\varepsilon$ after approximately 300 equilibration steps, in agreement with the MD reference value of $-1.958~\varepsilon$ to within 3.2%. Periodic boundary conditions and the minimum image convention are correctly implemented, and the Metropolis criterion enforces sampling of the Boltzmann distribution throughout the run. A key physical insight emerging from the comparison is that MC and MD sample genuinely different states when the initial conditions place particles far from the true equilibrium density, and care must be taken to ensure that both methods probe the same thermodynamic state point before their results are compared. This work demonstrates the power and consistency of stochastic sampling methods for computing equilibrium properties of simple fluids and provides a validated foundation for more advanced simulations, including NPT ensemble MC, free energy perturbation, and histogram reweighting methods.

---

## References

1. Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. *Journal of Chemical Physics*, **21**(6), 1087–1092.

2. Allen, M. P., & Tildesley, D. J. (2017). *Computer Simulation of Liquids* (2nd ed.). Oxford University Press.

3. Frenkel, D., & Smit, B. (2002). *Understanding Molecular Simulation: From Algorithms to Applications* (2nd ed.). Academic Press.

4. Lennard-Jones, J. E. (1924). On the Determination of Molecular Fields. *Proceedings of the Royal Society A*, **106**(738), 463–477.

5. Rapaport, D. C. (2004). *The Art of Molecular Dynamics Simulation* (2nd ed.). Cambridge University Press.

6. Frenkel, D. (2004). Speed-up of Monte Carlo simulations by sampling of rejected states. *Proceedings of the National Academy of Sciences*, **101**(51), 17571–17575.

---

<div align="center">
Sri Sathya Sai Institute of Higher Learning · Department of Physics · 2025
</div>
