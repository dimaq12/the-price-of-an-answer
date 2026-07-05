# The Price of an Answer

*Spectra Without Matrices · Journey II. Part I — [Never Quantum at All](https://github.com/dimaq12/do-we-need-quantum-computing) — used one matrix-free number to tell genuine quantum speedups from costumes. This time the same number walks into the world of nonlinear PDEs and asks every solver the rudest possible question: what am I actually paying you for?*

## The bill

Every answer in computational science arrives with a bill. You want to know how a shock forms at viscosity ν = 0.031? Run the solver. At 0.032? Run it again. A parameter sweep is a thousand runs; a calibration loop is a thousand sweeps. The solver charges per answer, and everyone pays, because the equations are *nonlinear* — and nonlinear means no shortcuts. Right?

I took 35 famous nonlinear PDEs — Burgers, KdV, Kuramoto–Sivashinsky in full chaos, a bright soliton, Camassa–Holm, fractional heat, thirty more — and put that assumption on a stand. The question, same as in Journey I: how much *structure* does the problem have, and what does that structure make redundant? There, the redundant machine was a quantum computer. Here it is the solver itself.

## One number, again

The dial is the effective rank Φ₁ — the number of directions an object really uses. In Journey I it read the spectra of data matrices. Here it reads the *parametric manifold*: solve the PDE at a handful of parameter values, stack the solutions, and ask the SVD how many independent directions the family {u(T; k)} actually spans.

For 33 of the 35 equations the answer is **Φ₁ = 1.00**. Chaos included. One direction. The entire family of answers, across a ±30% parameter range, is a one-dimensional curve bending through a 128-dimensional space.

Low rank meant "dequantizable" in Journey I. Here it means something just as blunt: *the solver is redundant at query time*. You do not need to re-solve a nonlinear PDE to move along a one-dimensional curve. You need to interpolate.

## Seven solves buy everything

The stand solves each equation at 7 Chebyshev points of the parameter range and builds a barycentric interpolant — a fifty-line construction from a numerical analysis textbook. That is the whole method. The price list it produces:

- **Forward answers:** 34 of 35 equations at relative error ~5·10⁻¹⁵ — machine precision — *across the entire parameter range*, not near an anchor. Queries take 8–14 microseconds: a median **576× faster** than the solver, and the ratio is honest — it includes nothing amortized, no GPU, no training.
- **Budget:** median 10 solves per equation (7 nodes + 3 held-out validation solves). The stand escalates only where the parametric curve is genuinely stiff — fractional heat, with the parameter in the *exponent* of the operator, needed 31 nodes. The stand reports every solve it spends.
- **Full disclosure:** we first built this kernel from Taylor sensitivities (W = ∂u/∂k — the tangent linear model), then tried to kill it with the strongest classical baseline we could think of. The baseline — Chebyshev interpolation on the same budget — won by five orders of magnitude. So we absorbed it. What you are reading is the method that beat our method.

## The same seven solves, read backwards

An interpolant you can query in microseconds is an inverse solver you can run for free. Hide a parameter k*, hand me one noisy snapshot of the solution, and finding k̂ is a scan along the one-dimensional curve — **zero additional PDE solves**.

- Clean data: **35 of 35** recovered below 10⁻¹⁰; several to the exact float64 bit.
- Noisy data: the estimator's median error lands at ~0.7× the Cramér–Rao bound — exactly where the theoretical 0.674 of an *optimal* estimator sits. This is the strongest sentence in the piece and it is worth saying plainly: past this point, better methods do not exist; the remaining error belongs to the noise, not the algorithm.
- The derivative of the interpolant is the Fisher information, so every estimate ships with an honest error bar — computed from the data itself. When a parameter is unidentifiable, the bar says so instead of the estimate lying.

## A whole map from one snapshot

The rudest version of the inverse problem: hide not a number but a *function* — a viscosity map ν(x), background plus a smooth wave plus a local defect — and hand over a single snapshot of Burgers flow at 1% noise. 128 unknowns, 128 noisy observations, and Hadamard laughing in the background.

The dial goes first: Φ₁ of the sensitivity operator reads **5.2**, and the singular spectrum says exactly **10 of 128** directions rise above the noise. One snapshot contains about ten independent numbers about the map. So the stand reconstructs precisely those ten — finds the defect at x = 4.27 (truth: 4.20) — and for everything else draws a wide band that says *the data are silent here*. It does not hallucinate the invisible 70%. And the bands are honest: 99.2% of the true map lies inside them.

Two lines from the stand deserve framing. Raise the noise to 5% and the visible count drops to **zero** — that snapshot knows nothing about the map, and the method says so out loud. Add a second snapshot and the count climbs to 12, the error falls — you can *buy* visibility with data, and the dial quotes the exchange rate.

This is the part that classical field inversion and neural reconstructions both get wrong by design: they always return a map, confident everywhere, half of it invented. Measuring what the data can see *before* reconstructing is the entire difference between an answer and a guess.

## The walls — and the dial saw them first

Honesty section, as in Journey I. One thing refused to fall — and the refusal was *predicted*.

The bright NLS soliton is the one equation whose forward interpolation resists: we escalated to 63 Chebyshev nodes — 200 solves — and the error crawled from 2·10⁻² only to 7·10⁻³. Its parametric family genuinely bends — the profile breathes, and moving structures do not compress into linear bases (the Kolmogorov n-width barrier, a theorem-grade wall). The dial flagged it before we knew: **Φ₁ = 2.70**, the only reading far from 1.00 in the whole set. The soliton is this journey's Shor: same dial, same verdict — *no cheap handle here, the structure is real*. (Its parameter, for the record, still inverts to the exact bit — through the amplitude, the one observable that stays tame.)

One more thing deserves to be said out loud, because it turns the title from a slogan into a theorem's shadow. A century ago Carleman showed that any polynomial nonlinearity becomes *linear* in a lifted basis of monomials — exactly linear, no approximation, just more coordinates (resona ships the construction: `resona.lift.carleman`). Every rank-1 parametric response in the table above is the visible edge of such a lift: a finite chart in which the nonlinear family straightens out. The soliton is the one case whose chart refuses to stay finite — and that refusal *is* the Kolmogorov n-width wall the dial flagged. Φ₁ has measured the same thing all along, in both journeys: the size of the linear chart you would need.

Scope, stated plainly: the kernel is built per initial condition and per time horizon; the test horizons are short; the noiseless inverse targets were generated by the same discretization (no model error); the components — Chebyshev interpolation, sensitivity analysis, maximum likelihood — are all classical. We did not invent new mathematics. We measured where the expensive machine is redundant, and it turned out to be almost everywhere we looked.

And the neural-operator elephant: the same parametric task is a standard benchmark for operator-learning papers, which typically spend thousands of training solves to reach 10⁻²–10⁻³. Our stand spends ten solves and reaches 10⁻¹⁴, plus error bars, plus a diagnosis. We leave that duel for the reader to run — the stands print every number.

## So — what does an answer cost?

Not what the solver charges. The honest price of an answer is set by the structure of the question: Φ₁ ≈ 1 and machine precision costs seven solves and a microsecond; Φ₁ = 2.7 and you have met a genuine wall; Φ₁ of your *data* is 10, and no method on earth will extract an eleventh number from them. Measure first. Then pay.

In Journey I the measured verdict was that some celebrated speedups were never quantum at all. This journey ends the same way, one shelf lower: within the range your questions live in, almost every famous nonlinear PDE was — *never nonlinear at all*.

## Run it yourself

Everything — this essay, all three stands, and the figures — lives in one repository: **[github.com/dimaq12/the-price-of-an-answer](https://github.com/dimaq12/the-price-of-an-answer)**.

```
git clone https://github.com/dimaq12/the-price-of-an-answer.git
cd the-price-of-an-answer
pip install -r requirements.txt
python stands/seven_solves.py       # the full price list over 35 PDEs
python stands/crb_optimal.py        # Cramér–Rao-optimal inversion, error bars included
python stands/one_snapshot_map.py   # the map from one snapshot, blind zones measured
```

Every number above is printed by a stand, checked against a hidden truth. And the machinery is now a library: **pip install lastsolve** — the accelerator, the dial, the certificates and the honest refusals, packaged ([github.com/dimaq12/lastsolve](https://github.com/dimaq12/lastsolve)). Don't trust me — run it.
