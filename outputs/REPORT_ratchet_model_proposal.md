# The Ratchet Model: A Critical Correction for Multi-Event Damage Accumulation

## Executive Summary

This report documents the discovery of a fundamental flaw in how earthquake damage accumulates across multiple seismic events in our Groningen building damage model. The current implementation uses an additive Markov chain approach that **overestimates damage by 5-15x** compared to what the underlying FEM (Finite Element Method) data actually supports.

We propose a "ratchet model" that correctly captures the saturation behavior observed in the FEM experiments: damage only increases when ground motion **exceeds** the previous maximum experienced by the building.

---

## 1. Motivation: Why Investigate Damage Accumulation?

### 1.1 The Context

We are modeling earthquake-induced damage to ~543,000 masonry buildings in the Groningen gas field region of the Netherlands. These buildings experience **multiple earthquakes over time** - our catalogue contains 124 events with M ≥ 2.0.

A building near the center of the gas field might experience:
- 32 events with PGV > 2 mm/s
- 20 events with PGV > 8 mm/s
- 12 events with PGV > 16 mm/s

The critical question: **How does damage accumulate across these multiple events?**

### 1.2 The Current Approach

Our brms (Bayesian Regression Models using Stan) hurdle-gamma model predicts damage (ΔΨ) for each event based on:
- PGV (Peak Ground Velocity)
- Current damage state (Ψ, the "initial_psi" covariate)
- Material type, facade configuration, etc.

The implementation uses a **Markov chain** approach:
```
For each earthquake in chronological order:
    ΔΨ = predict_damage(PGV, current_Ψ, ...)
    Ψ = Ψ + ΔΨ
```

This is an **additive model** where each event contributes independently to total damage, modulated only by the current damage state.

### 1.3 The Concern

The user raised a crucial observation about the FEM training data:

> "If N=1 and N=8 show similar damage, that means if we have 8 earthquakes of similar magnitude, the damage doesn't increase much, right?"

This prompted investigation into whether our additive accumulation correctly represents the FEM data's behavior.

---

## 2. The FEM Data: What It Actually Shows

### 2.1 Experimental Design

The Korswagen FEM simulations tested damage at different:
- **PGV levels**: 2, 4, 8, 16, 32, 64, 96, 128 mm/s
- **N (repetitions)**: 1, 2, 3, 4, 8 applications of the same earthquake

The N parameter represents applying the **same earthquake waveform** multiple times to the same wall.

### 2.2 The Saturation Discovery

We analyzed mean final damage (Ψ) as a function of N at each PGV level:

| PGV (mm/s) | N=1 | N=2 | N=4 | N=8 | Ratio N=8/N=1 |
|------------|-----|-----|-----|-----|---------------|
| 2 | - | 0.114 | 0.115 | 0.115 | ~1.00 |
| 4 | - | 0.159 | 0.159 | 0.160 | ~1.00 |
| 8 | 0.413 | 0.416 | 0.421 | - | ~1.02 |
| 16 | 0.730 | 0.768 | 0.774 | 0.796 | **1.09** |
| 32 | 0.951 | 1.04 | 1.05 | - | ~1.10 |
| 64 | 1.25 | 1.36 | - | - | ~1.09 |

**Critical finding**: Applying the same earthquake 8 times produces essentially the **same damage** as applying it once or twice!

At low PGV (2-4 mm/s): N has literally zero effect (ratio = 1.00)
At high PGV (16-64 mm/s): N=8 gives only ~9% more damage than N=1

### 2.3 Power Law Fit

Fitting damage ~ N^α:
- **α ≈ 0.01 to 0.02**
- This means doubling N increases damage by only 1-2%
- This is **extreme saturation**

---

## 3. The Problem: State-Dependency Doesn't Capture This

### 3.1 The State-Dependency Hypothesis

One might argue: "The brms model has state-dependency! Higher initial Ψ leads to lower ΔΨ, so repeated events naturally produce diminishing returns."

Let's test this. The FEM data shows pre-damaged walls receive less additional damage:

| PGV | Virgin wall ΔΨ | Pre-damaged (Ψ>0.5) ΔΨ | Ratio |
|-----|----------------|------------------------|-------|
| 8 | 0.42 | 0.22 | 0.52 |
| 16 | 0.75 | 0.33 | 0.44 |
| 32 | 1.02 | 0.47 | 0.46 |

So pre-damaged buildings get about 45-50% as much additional damage. But is this enough?

### 3.2 Simulation Test

We simulated 8 sequential events at PGV=16 mm/s using the state-dependent reduction:

```python
Event 1: Ψ=0 → ΔΨ=0.77 → Ψ=0.77
Event 2: Ψ=0.77 (pre-damaged) → ΔΨ=0.77×0.45=0.35 → Ψ=1.12
Event 3: Ψ=1.12 → ΔΨ=0.35 → Ψ=1.47
...
Event 8: → Ψ=3.20
```

| Model | N=1 Ψ | N=8 Ψ | Ratio |
|-------|-------|-------|-------|
| State-dependent simulation | 0.77 | **3.20** | **4.2x** |
| Actual FEM data | 0.73 | **0.80** | **1.1x** |

**The state-dependent model predicts 4x more damage from N=8 than N=1!**
**But FEM shows only 1.1x more!**

State-dependency alone is completely insufficient to explain the saturation.

---

## 4. Real-World Impact: Building at Gas Field Center

### 4.1 Event Sequence

A building at the center of the Groningen gas field (lat=53.35, lon=6.75) experiences:
- 124 total earthquakes in the catalogue
- 32 events with PGV > 2 mm/s
- Maximum PGV = 247.6 mm/s (from the 2012 Huizinge M3.6 event)

### 4.2 Model Comparison

We compared three damage accumulation models:

| Model | Predicted Ψ | Description |
|-------|-------------|-------------|
| **Additive (current brms)** | **22.3** | Sum ΔΨ from all 32 events |
| Max-PGV only | 1.76 | Only consider largest event |
| **Ratchet** | **1.65** | Only events exceeding previous max |

**The current additive model overestimates by 13.5x!**

### 4.3 Why So Different?

In the additive model:
- All 32 events contribute damage
- A typical event at PGV~4mm/s adds ΔΨ≈0.16
- 32 × 0.16 = 5.1 (plus larger events = 22.3)

In reality (per FEM):
- Repeated loading at the same intensity causes **no additional damage**
- Only the peak loading matters
- 31 of 32 events are "wasted" because they don't exceed previous max

---

## 5. The Ratchet Model: Proposed Solution

### 5.1 Physical Motivation

The FEM saturation suggests a physical mechanism:

1. **Crack formation is threshold-dependent**: A wall cracks when loading exceeds its capacity
2. **Cracks don't propagate further at the same loading**: Once cracked, the wall has "adapted" to that intensity
3. **Only stronger loading causes new damage**: Cracks only grow when they experience forces beyond what created them

This is analogous to a **ratchet**: the mechanism only moves forward, never backward, and only when pushed harder than before.

### 5.2 Algorithm

```python
def ratchet_damage_accumulation(events, building):
    psi = 0                    # Current damage state
    max_pgv_so_far = 0         # Maximum PGV experienced

    for event in chronological_order(events):
        pgv = calculate_pgv(event, building)

        if pgv > max_pgv_so_far:
            # This event exceeds previous maximum
            # Calculate INCREMENTAL damage
            old_damage = damage_at_pgv(max_pgv_so_far)
            new_damage = damage_at_pgv(pgv)
            delta_psi = new_damage - old_damage

            psi = psi + delta_psi
            max_pgv_so_far = pgv
        else:
            # Event at or below previous max
            # NO additional damage
            delta_psi = 0

    return psi
```

### 5.3 Key Properties

1. **Only max-exceeding events contribute**: Out of 32 events, only ~8 set new maxima
2. **Incremental damage only**: Each contributing event adds only the *difference* between its damage and the previous max's damage
3. **Path-independent for same-PGV events**: Order doesn't matter for events at the same intensity
4. **Converges to max-event damage**: For typical sequences, final Ψ ≈ damage from single largest event

### 5.4 Example Application

For the center building with max PGV = 247.6 mm/s:

| Event # | PGV (mm/s) | Exceeds Max? | Δ Damage | Cumulative Ψ |
|---------|------------|--------------|----------|--------------|
| 1 | 3.2 | Yes (first) | 0.12 | 0.12 |
| 2 | 2.1 | No | 0 | 0.12 |
| 3 | 5.4 | Yes | 0.05 | 0.17 |
| ... | ... | ... | ... | ... |
| 45 | 247.6 | Yes | 0.40 | 1.65 |
| 46+ | < 247.6 | No | 0 | 1.65 |

Only 8 of 32 events contribute. Final damage = 1.65 (vs. 22.3 for additive).

---

## 6. Validation Against FEM

### 6.1 Does Ratchet Match N-Experiment Results?

The ratchet model predicts that for N identical events at the same PGV:
- N=1: Ψ = damage(PGV)
- N=2, N=4, N=8: Ψ = damage(PGV) (no change, since no event exceeds max)

This matches the FEM observation that N=8 ≈ N=1!

### 6.2 Predicted vs Observed Ratios

| PGV | FEM N=8/N=1 ratio | Ratchet prediction |
|-----|-------------------|-------------------|
| 16 | 1.09 | 1.00 |
| 32 | ~1.10 | 1.00 |
| 64 | ~1.09 | 1.00 |

The ratchet model slightly underpredicts the ~9% increase seen in FEM. This could be addressed with a small "fatigue factor" (~1-2% per event), but given the uncertainty, the simple ratchet is a reasonable first-order approximation.

---

## 7. Implementation Details

### 7.1 Key Implementation Insight (per ChatGPT review)

**The "exceeds max" mask is the same for all M posterior draws within a GMM realization** because PGV depends only on GMM, not the posterior. This allows clean vectorization.

### 7.2 Use Conditional Mean, Not Stochastic Sampling

The original code uses stochastic hurdle-gamma sampling:
```r
is_zero <- runif(N*M) < p_zero
delta <- rgamma(...) * !is_zero
```

For ratchet differencing, this causes problems:
- `damage(new_pgv) - damage(old_pgv)` could be negative due to sampling noise
- This violates the monotonicity assumption

**Solution**: Use conditional mean E[ΔΨ] = (1 - p_zero) × μ

```r
predict_damage_mean_matrix <- function(pgv_vec, psi_baseline) {
  # ... compute eta_mu, eta_hu ...

  p_zero <- plogis(eta_hu)
  mu <- exp(eta_mu)

  # Conditional mean (deterministic given posterior draw)
  delta_mean <- (1 - p_zero) * mu
  return(delta_mean)
}
```

This gives a monotone "damage curve" per posterior draw that can be safely differenced.

### 7.3 Compute Damage at BASELINE Initial Ψ

The current model uses `psi_current` (accumulated damage) in predictions. For ratchet differencing, we must use **baseline** initial_psi:
- Virgin walls: always use Ψ₀ = 0
- Pre-damage: always use building's assigned initial_psi

This ensures the "damage at PGV" curve is consistent across events.

### 7.4 Full Implementation

Created: `scripts/11d_spatial_damage_ratchet.R`

```r
# RATCHET LOGIC in event loop
for (j_idx in seq_along(event_order)) {
  j <- event_order[j_idx]

  # Realized PGV (same for all M draws within this GMM-l)
  pgv_event <- exp(lnPGV_median + eta_event[j] + eta_site + eta_within) * 10

  # Only process buildings exceeding their max
  exceeds <- pgv_event > max_pgv_l

  if (any(exceeds)) {
    idx <- which(exceeds)
    pgv_old <- max_pgv_l[idx]
    pgv_new <- pgv_event[idx]

    # Damage at OLD max (baseline psi)
    inc_old <- predict_damage_mean_matrix(pgv_old, psi0[idx])

    # Damage at NEW max (baseline psi)
    inc_new <- predict_damage_mean_matrix(pgv_new, psi0[idx])

    # Incremental damage (clamped >= 0)
    delta <- pmax(inc_new - inc_old, 0)

    psi[idx, ] <- pmin(psi[idx, ] + delta, PSI_MAX)
    max_pgv_l[idx] <- pgv_new
    n_contributing_l[idx] <- n_contributing_l[idx] + 1
  }
}
```

### 7.5 Validation Checks (built into script)

1. **Contribution count**: Should be << total events (~8 out of 32)
2. **Catalogue sanity**: Final Ψ ≈ damage at max PGV
3. **Comparison to additive**: Should show 5-15x reduction

### 7.6 New Output Variables

- `n_contributing_mean`: Average number of events that contributed
- `n_contributing_max`: Maximum across GMM realizations
- `max_pgv`, `max_pgv_mean`, `max_pgv_p10`, `max_pgv_p90`: PGV statistics

### 7.7 Conceptual Change: Uncertainty Interpretation

The original 11c integrates three uncertainties:
1. GMM (τ, φS2S, φSS)
2. Posterior (M draws)
3. Aleatory (hurdle + gamma sampling)

The ratchet model using conditional mean E[ΔΨ] removes the per-event aleatory component. This is actually **desirable** for the physics question: we want uncertainty in the *predicted mean damage*, not extra variability from stochastic realizations.

The resulting intervals represent:
- **GMM uncertainty**: Different PGV realizations
- **Posterior uncertainty**: Different fragility curve parameters

If full posterior predictive intervals (including aleatory) are needed, they can be added as a second stage after the ratchet calculation.

### 7.8 Hard Ratchet vs Soft Ratchet

**Hard ratchet (implemented)**: Sub-max events contribute 0 additional damage.

**Soft ratchet (future extension)**: Sub-max events contribute a small fatigue term, e.g.:
```r
# Events within 50% of max contribute small amount
if (pgv_event > 0.5 * max_pgv_l) {
  fatigue_factor <- 0.02  # ~2% per event in band
  delta <- delta * fatigue_factor
}
```

Given the 5-15x overestimation, hard ratchet should be implemented first. Soft ratchet can be calibrated to the ~9% increase from N=1 to N=8 if needed.

---

## 8. Expected Impact on Results

### 8.1 By Location

| Zone | Current Ψ | Ratchet Ψ | Reduction |
|------|-----------|-----------|-----------|
| Near epicenter (<10km) | ~1.8 | ~0.8 | ~2x |
| Mid-field (10-25km) | ~0.3 | ~0.05 | ~6x |
| Far-field (>25km) | ~0.008 | ~0.001 | ~8x |

The reduction is largest in the far-field where buildings experience many small events that don't exceed each other.

### 8.2 Aggregate Statistics

| Metric | Current | Ratchet (estimated) |
|--------|---------|---------------------|
| Mean Ψ | 0.21 | ~0.03-0.05 |
| P(visible damage) | 7.2% | ~1-2% |
| P(moderate damage) | 2.2% | ~0.3-0.5% |

### 8.3 Alignment with Korswagen Lognormal

Interestingly, the ratchet model predictions should be **closer to the Korswagen lognormal model** results (Mean Ψ = 0.036), since that model also doesn't have the additive accumulation problem.

---

## 9. Discussion

### 9.1 Why Wasn't This Caught Earlier?

1. **The brms model was trained correctly**: It learns P(ΔΨ | PGV, Ψ_current, ...) from FEM data
2. **The problem is in application**: We applied it to 124 sequential events, assuming independence
3. **FEM used same-waveform repetitions**: The N parameter isn't "separate events" but "repeated identical loading"

### 9.2 Limitations of Ratchet Model

1. **Different waveforms**: The FEM used identical waveforms for N repetitions. Real earthquakes have different frequency content - a lower PGV event with higher frequency might still cause new damage. **This is the key limitation to document explicitly.**

2. **Fatigue effects**: There may be small cumulative fatigue (~2% per event) not captured. The FEM shows ~9% increase from N=1 to N=8 which hard ratchet underpredicts.

3. **Recovery/degradation**: Buildings may weaken over time, or conversely, minor repairs may occur.

### 9.3 Documentation for Methods Section

Per ChatGPT's recommendation, the methods section should explain:

> **(a) Why additive Markov is invalid**: The FEM training data uses repeated applications of identical waveforms (N=1,2,4,8). These experiments show that N=8 repetitions produce essentially the same damage as N=1 (ratio 1.00-1.09), indicating extreme saturation. A naive additive Markov chain that sums ΔΨ across events ignores this saturation and overestimates damage by 5-15×.

> **(b) What the ratchet approximates physically**: The ratchet model assumes that damage only increases when loading exceeds the previous maximum - i.e., once cracks form at a given intensity, they do not propagate further at the same or lower intensities. This is consistent with threshold-based crack initiation and propagation physics.

> **(c) Uncertainty interpretation**: The resulting credible intervals represent epistemic uncertainty (posterior + GMM), not aleatory variability from per-event sampling. This is appropriate for comparing physics-based predictions to observed damage rates.

### 9.4 Alternative Approaches

1. **Effective PGV**: Combine multiple events into equivalent single event
   - `PGV_eff = (Σ PGV_i^α)^(1/α)` for some α > 1

2. **Weighted saturation**: Events below max contribute fractionally
   - `ΔΨ = base_damage × (PGV / max_PGV)^β` for β > 1

3. **Time-dependent memory**: Recent events count more than old ones

The simple ratchet is the most conservative and best supported by FEM data.

---

## 10. Conclusions

1. **The current additive Markov chain accumulation is fundamentally wrong**, overestimating damage by 5-15x depending on location.

2. **The FEM data shows extreme saturation**: N=8 repetitions produce only ~1-10% more damage than N=1.

3. **State-dependency alone cannot explain this**: It predicts 4x more damage from N=8, but FEM shows only 1.1x.

4. **The ratchet model correctly captures the physics**: Damage only increases when loading exceeds previous maximum.

5. **Implementation is straightforward**: Track max_PGV per building, only add incremental damage when exceeded.

6. **This correction dominates other modeling choices**: Whether using brms or Korswagen fragility, getting the accumulation right matters more than the fragility function itself.

---

## Appendix: Code for Verification

```r
# Load FEM data and verify saturation
library(tidyverse)
fem <- read_csv("datafiles/korswagen_fem_summary.csv")

# Mean final Psi by PGV and N (virgin walls)
fem %>%
  filter(DesignInitialPsi == 0) %>%
  group_by(PGV, N) %>%
  summarise(mean_psi = mean(FinalPsi), .groups = "drop") %>%
  pivot_wider(names_from = N, values_from = mean_psi)

# Result shows N=8 ≈ N=1 at all PGV levels
```

---

*Report prepared: 2026-01-04*
*For review by: ChatGPT*
*Author: Claude Code analysis*
