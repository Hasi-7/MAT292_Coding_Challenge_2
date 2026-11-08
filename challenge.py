"""Clock recovery: edit a parameter vector and compare it with the target."""
import numpy as np
from challenge_utils import (
    SCORE_IMPROVEMENT_THRESHOLD,
    evaluate_parameters,  # (parameters) -> objective and components
    plot_model,           # (parameters, filename="challenge_clock.png") -> saves a plot
    print_objective,      # (result, parameters, title) -> prints a formatted report
    simulate,
)
from constants import REPORTER_PARAMETERS, GROUP_INDICES
from scipy.optimize import differential_evolution

def cost(q, p_rep, indices):
    p_trial = p_rep.copy()
    p_trial[indices] = q
    return evaluate_parameters(p_trial)["objective"]

def candidate(genes, p_rep = REPORTER_PARAMETERS.copy()):
    indices = gene_indices(genes)
    bounds = [
    (0.01 * REPORTER_PARAMETERS[i], 2.5 * REPORTER_PARAMETERS[i])
    for i in indices
    ]

    result = differential_evolution(
        cost, 
        bounds=bounds,
        args=(p_rep, indices,),
        rng=np.random.default_rng(123))

    if not result.success:
        raise RuntimeError(result.message)
    
    p_edit = p_rep.copy()
    p_edit[indices] = result.x

    return p_edit


# Functions imported from challenge_utils:
#
# evaluate_parameters(parameters)
#     Simulate the clock for a 45-entry parameter vector and return a dict
#     with the objective J and its components, using the supplied period estimator.
#
# plot_model(parameters, filename="challenge_clock.png")
#     Plot the PER, CRY1, CRY2, ROR and REV trajectories and save the
#     figure to filename.
#
# print_objective(result, parameters=REPORTER_PARAMETERS, title="Objective report")
#     Print a formatted report. result is the output of evaluate_parameters,
#     parameters is the vector that produced it, and title is the heading.


def gene_indices(genes):
    """Return full-vector indices belonging to the named genes, in gene order."""
    if isinstance(genes, str):
        genes = [genes]
    return np.concatenate([GROUP_INDICES[gene] for gene in genes])


def update_pars(pold, pnew, genes):
    """Copy selected gene entries from reduced pnew into full vector pold."""
    updated = np.asarray(pold, dtype=float).copy()
    indices = gene_indices(genes)
    if np.asarray(pnew).shape != (len(indices),):
        raise ValueError(f"pnew must contain {len(indices)} values for {genes}.")
    updated[indices] = pnew
    return updated

if __name__ == "__main__":
    # Keep the reporter reference unchanged.
    p_rep = REPORTER_PARAMETERS.copy()

    # This is the part where you would create your own p_edit
    # An outline of how to use gradient descent on a selected gene set is given below
    # The gradient_calculator is the function you would need to write
    # using the ideas from Q3.
    GENES = ["CRY1"]
    p_edit = candidate(GENES)

    result = evaluate_parameters(p_edit) 
    print_objective(result, p_edit, "OPTIMIZED CANDIDATE") # Score the final vector after your algorithm.
    plot_model(p_edit)
    print(p_edit)
    # Check the final score against the assignment threshold (higher is better). 
    score_improvement = result["score_improvement_percent"]  # Percent improvement in the full score J.
    if score_improvement is None:
        print("The final vector could not be scored: " + result.get("failure", "simulation failed"))
    elif result.get("valid", False) and score_improvement >= SCORE_IMPROVEMENT_THRESHOLD:
        print(f"You are above the assignment threshold ({score_improvement:.1f}% improvement).")
    else:
        print(f"You are below the assignment threshold ({score_improvement:.1f}% improvement).")
