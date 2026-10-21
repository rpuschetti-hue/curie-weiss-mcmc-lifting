# Non-Reversible MCMC: Lifted Metropolis-Hastings on the Curie-Weiss Model

This repository contains the MATLAB numerical implementation developed for my Master's Thesis in Applied Mathematics at **Sapienza University of Rome**: *"Non-Reversible Markov Chains: Quantitative Analysis and Convergence Acceleration"*.

## 📌 Overview & Mathematical Background
Standard reversible Markov Chain Monte Carlo (MCMC) algorithms, such as the classic Metropolis-Hastings (MH), suffer from severe critical slowing down near phase transitions due to their diffusive nature.

By applying **State-Space Lifting** (following the framework of *Turitsyn, Chertkov, Vucelja* and the scaling limits of *Bierkens & Roberts*), the state space $\Omega$ is doubled into $\Omega \times \{+1, -1\}$ to introduce directional persistence and break detailed balance. 

At the critical inverse temperature ($\beta = 1$), this non-reversible dynamics achieves a quadratic-like speedup in the relaxation time $\tau$:
* **Reversible Metropolis-Hastings:** $\tau \sim \Theta(N^{3/2})$ (Diffusive Langevin limit)
* **Lifted Metropolis-Hastings:** $\tau \sim \Theta(N^{3/4})$ (Ballistic Piecewise Deterministic Markov Process / Zig-Zag limit)

## 📊 What the Code Does
The script `curie_weiss_lifting_analysis.m` performs a complete spectral and stochastic analysis of the macroscopic magnetization chain ($N$ spins), generating **5 diagnostic figures**:

1. **Asymptotic Relaxation Scaling:** Log-log regression comparing singular value gaps across system sizes $N \in [20, 640]$, confirming the $\Theta(N^{1.5})$ vs. $\Theta(N^{0.75})$ theoretical slopes.
2. **Temporal Error Decay:** Evolution of the squared expectation $\langle m(t) \rangle^2$ starting from a polarized state ($m_0 = +1$), showcasing diffusive decay versus non-reversible damped oscillations.
3. **Topology Schematic:** Visual representation of reversible random walks versus lifted persistent flows.
4. **Stochastic Trajectories (Zig-Zag Limit):** Comparison of sample paths at equilibrium, highlighting the emergence of Piecewise Deterministic Markov Process (PDMP) trajectories.
5. **Complex Eigenvalue Spectrum:** Visualization of the transition matrix eigenvalues in the complex unit circle, proving the breaking of self-adjointness and detailed balance.

## 🚀 How to Run
1. Clone this repository or download `curie_weiss_lifting_analysis.m`.
2. Open **MATLAB** (no external toolboxes required).
3. Run the script:
   ```matlab
   run('curie_weiss_lifting_analysis.m')
