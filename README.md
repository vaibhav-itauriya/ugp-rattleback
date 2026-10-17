# Slipping-contact friction in rattleback dynamics

**Extension of a rolling formulation using the Moffatt–Shimomura contact laws**

Vaibhav Itauriya and Dipayan Mukherjee
Department of Mechanical Engineering, Indian Institute of Technology Kanpur
Undergraduate project (UGP)

📄 **Paper:** [`paper/Itauriya_Mukherjee_rattleback.pdf`](paper/Itauriya_Mukherjee_rattleback.pdf)

🌀 **Interactive simulator:** <https://vaibhav-itauriya.github.io/ugp-rattleback/>. Choose the contact model, friction and launch conditions; it computes the motion in your browser and shows the results, all the time histories and a 3D view of the rattleback.

---

## Summary

A rattleback spun in its non-preferred sense converts spin into a rocking
oscillation and comes out spinning the other way. Most models of this assume
rolling without slipping, which removes friction and leaves the contact force
unresolved. This project releases the rolling constraint (6 → 9 first-order
equations), builds the contact force explicitly, and solves the normal reaction
in closed form from the condition that the body stays on the plane:

$$
N=\frac{m\,(g+c+\mathbf d\cdot\mathbf A_0)}{1-m\,\mathbf d\cdot\mathbf B-m\,\mathbf d\cdot\mathbf A_1}
$$

Three tangential laws are compared: regularised Coulomb, viscous (the two laws of
Moffatt & Shimomura, *Nature* **416**, 385, 2002), and Coulomb with a drilling torque.

| Contact law (ω₀ = −5 rad/s, 8 s) | Reversal | θ<sub>max</sub> | ω₃(8 s) | ΔE/E₀ |
|---|---|---|---|---|
| Rolling | 2.83 s | 6.48° | +5.00 | 0 |
| Coulomb, μ = 0.15 | 4.22 s | 0.84° | +1.01 | −12.4% |
| Coulomb + drilling, ρ = 0.5 mm | 3.98 s | 0.84° | +0.02 | −12.9% |
| Viscous, c = 500 kg/s | 2.90 s | 3.09° | +1.40 | −11.9% |

Main findings:

- In the rolling solution the normal reaction ranges from 0.20 mg to 5.86 mg.
  Holding it at mg shifts the reversal time by 14% near the rolling limit and
  prevents convergence to that limit.
- A viscous law of the same nominal strength as Coulomb friction produces no
  reversal at all, mirroring the finite-time versus asymptotic transitions of
  Moffatt & Shimomura.
- A drilling torque with a Hertzian lever arm (ρ ≈ 0.16–0.52 mm) keeps the
  reversal and then spins the body down; at ρ = 2 mm it stops the body before
  the instability can grow.

## Repository layout

```
paper/
  Itauriya_Mukherjee_rattleback.pdf   the paper
  rattleback.tex, figs/               LaTeX source (revtex4-2) and figures
code/
  matlab/
    rattleback_simulator_123H.m       main code: all models, sweeps, figures
    console_output.txt                output of the run used for the paper
    figures/                          figures exported by that run
  python/                             independent cross-check and movie renderer
movies/                               six animations (0.4x real time)
docs/                                 interactive simulator (GitHub Pages site)
  js/model.js                         the same equations in JavaScript
tools/validate_model.mjs              checks model.js against the paper's numbers
```

## Running

**MATLAB** (R2021a or later):

```matlab
>> rattleback_simulator_123H
```

Figures are written to a `figures/` folder next to the `.m` file. Studies
(main, ε sweep, ρ sweep, μ sweep, initial spin, fixed-N comparison, GIF
animation) are switched on and off in `config()` at the top of the file. A full
run takes roughly 10–20 minutes.

**Python** (numpy, scipy, matplotlib, pillow; ffmpeg for movies):

```bash
cd code/python
python3 run_simulations.py          # runs -> runs/*.npz
python3 figures.py                  # schematic and cross-check figures
python3 videos.py coulomb_nonpref   # one movie; also rolling_nonpref,
                                    # coulomb_spin_nonpref, viscous_nonpref, grid, senses
```

The two implementations agree on reversal times to within 0.002 s.

**Browser simulator**: open the link above, or serve `docs/` locally
(`python3 -m http.server -d docs`) and visit <http://localhost:8000>. To check
its physics against the paper: `node tools/validate_model.mjs`.

## Parameters

m = 1 kg; a₁, a₂, a₃ = 0.20, 0.03, 0.02 m; b = 0.005 m;
I₁₁, I₂₂, I₃₃ = 2.0, 16.0, 17.5 × 10⁻⁴ kg m²; I₁₂ = −0.5 × 10⁻⁴ kg m²
(Kane & Levinson, *Int. J. Non-Linear Mech.* **17**, 175, 1982).
Regularisation ε = 3 × 10⁻⁵ m/s (converged; see the ε sweep).
