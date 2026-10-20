# Six completed representations under common downstream inference

The [declared exploratory plan](../config/norman_shared_inference_plan.json) was
fixed after original model/marker/phase outcomes and before the common-inference
outcomes. All six arms completed the declared protocol. The result preserves the
narrow conclusion: no SAE advantage or new biological regulator is established.

Every arm used the same 5,926 discovery and 5,929 inference controls. Discovery
controls alone supplied annotation; perturbation tests used inference controls
within gemgroup. All 6 × 9 × 105 = 5,670 association hypotheses entered joint BH,
including unmatched p=1 records. Calls require q<.05 and AUC≥.55. All 54 arm-panel
identity hypotheses entered joint Holm with 9,999 permutations. Each matched
panel had raw p=.0001 and joint Holm p=.0054; unmatched panels retain p=1.
Identity heterogeneity is not pathway specificity or independent replication.

| Panel | All-cell SAE 0 | SAE 1 | SAE 2 | Control-only SAE 0 | PCA | Capped NMF |
|---|---:|---:|---:|---:|---:|---:|
| Ery_TF | 0 | 0 | 5 | Unscored | 28 | 6 |
| Ery_Heme | Unscored | Unscored | Unscored | Unscored | Unscored | 21 |
| Ery_Membrane | Unscored | Unscored | Unscored | Unscored | Unscored | Unscored |
| Ery_Effector | 1 | 1 | 1 | 0 | Unscored | 54 |
| Gran_TF | Unscored | Unscored | Unscored | Unscored | Unscored | Unscored |
| Gran_Primary | Unscored | Unscored | Unscored | Unscored | Unscored | Unscored |
| Gran_Secondary | Unscored | Unscored | Unscored | Unscored | Unscored | Unscored |
| Progenitor | Unscored | Unscored | Unscored | Unscored | Unscored | Unscored |
| Cycling | 0 | 0 | 0 | 0 | 52 | 0 |

Numbers count joint association calls, not identified biological mechanisms.
The [signed-weight follow-up](NORMAN_PCA_ATTRIBUTION.md) explains the PCA mean
response without equating the score with a biological mechanism.
Unscored means no annotation met the fixed gate, not absence of biological effects.
PCA/NMF rank64 and SAE latent2048 have unequal capacity; NMF did not converge
within its original cap. Presence-based factor annotation also differs from direct
marker detection eligibility. Response counts cannot rank biological fidelity.

KLF1 is again the sole candidate appearing in at least two all-cell SAE seeds,
on the effector panel. This is computational agreement among models of the same
cells, and it does not repair the legacy background-gate failure. The seed-2 TF
list shrinks from 11 original legacy candidates to five under the new procedure;
no TF candidate repeats. The control-only heme annotation no longer passes
(adjusted match p=.0990), illustrating sensitivity to annotation controls rather
than a negative result for erythroid biology.

All six models match a Cycling feature. Only PCA has upward associations at the
fixed effect floor. Of its 52 calls, 32 have S AUC≤.45, 22 G2M AUC≤.45 and 20 both.
These are descriptive paired AUC directions without a new significance test.
The [separate phase diagnostic](NORMAN_PHASE_DIAGNOSTIC.md) reports its stricter
directional Holm results: among 50 retained PCA calls, 19 also have significant
decreases in both phases. The different counts answer different questions.

The original PCA fit_transform scores were reconstructed using the exact original
solver/settings and one-thread execution. Every one of its 945 archived AUCs,
p-values, match choices and control split reproduced within the declared strict
tolerances before their use. Actual saved NMF factors likewise reproduced their
archived inference. No SAE/NMF optimization was repeated. This avoids presuming
that a later linear PCA transform equals randomized-SVD fit_transform.

Complete [5,670 comparisons](../data_norman_continuation/shared_inference/comparisons.csv),
[summary](../data_norman_continuation/shared_inference/summary.json), complete
losslessly compressed JSON, executed log and
[validation record](../data_norman_continuation/shared_inference/validation.json)
are archived. Sources and actual activation hashes were checked. The all-cell
fits saw inference-cell expression, so the common split is conditional on fitted
representations, not a model-training holdout. The single control-only seed also
trained on inference controls. Shared cells, technical gemgroups and model seeds
are not independent biological samples; unequal severity and cell dependence
continue to limit inference.
