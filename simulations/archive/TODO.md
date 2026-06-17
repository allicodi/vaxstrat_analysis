- Add monotonicity sensitivity analysis to supplement
	- to supplement L?

- Update `do_gcomp_nat_inf` and `do_aipw_nat_inf` to use the correct formula for mu_1dot when both assumptions hold

- Simulations
	- "Asymptotic properties of estimators simulation" definitely needs to be re-run for all scenarios
		- There is actually an issue with the unadusted bounds calculation when P(Z | X) depends on X. I have submitted a PR where we now make `simulate_data_generic` look like a simple randomized trial.
		- So all tables for this sim need to be updated. 
		- Update text as needed. Mark changes. 
	- "Comparing power of estimands in realistic setting" simulation technically doesn't have to be re-run since both ER and PI assumptions hold. Thus, the semiparametric estimator is technically "correct" so I don't expect results will change much. But if not too much trouble, probably worth re-running in case we get slight improvement.
		- If you re-run, update text as needed if results change.
	- Need a new simulation to look at cross-fitted estimators. 
		- I submitted a PR with a new function that simulates data under a different DGP with various step functions and weird interactions. The structure is otherwise the same as previous functions.
		- I assume from a human-time-standpoint, it is easier to just keep the structure exactly the same as the generic simulation and go ahead and do it for all the estimators.
		- ONLY consider the PI and ER satisfied scenario. No need to violate those things. Set `doomed_inflation = 0.1` (or update text in Supp so that it's 0).
		- For the ML, you could just do random forest. Or if everything runs through super learner, then do SL.glm and SL.ranger or something. 
		- Just set cross-fitting folds to 3 or 5 or something small-ish so it runs quickly
		- Sample sizes 500 and 4000 are probably fine to start with
		- Add results to Supplement I.3. I have added some text there describing the DGP and left placeholders for true values (don't need a table since just one scenario)
		- Random questions: 
			- is [this line](https://github.com/allicodi/vaxstrat_analysis/blob/1749f15db20555ef023f925892bd5b4416f36055/simulations/R/run_simulation_1.R#L74) needed? Aren't we just reporting closed form inference?
			- `simulate_data_contour` is deprecated right? should we remove from repo?
		- We want to compare a cross-fitted estimator vs. non-cross fitted. So the output is like Table S3, except formatted as:

| assumptions satisfied | cross_fit | sample_size | bias | etc... |
|----------------------|----------|------------:|-----:|--------|
| PI                   | yes      | 400         | .... | ....   |
| PI                   | no       | 400         | .... | ....   |
| ER                   | yes      | 400         | .... | ....   |
| ER                   | no       | 400         | .... | ....   |
| both                 | yes      | 400         | .... | ....   |
| both                 | no       | 400         | .... | ....   |
| PI                   | yes      | 5000        | .... | ....   |
| PI                   | no       | 5000        | .... | ....   |
| ER                   | yes      | 5000        | .... | ....   |
| ER                   | no       | 5000        | .... | ....   |
| both                 | yes      | 5000        | .... | ....   |
| both                 | no       | 5000        | .... | ....   |

		
- Real data analysis
	- Covariate-adjusted bounds need to be re-run and semiparametric estimator needs to be re-run.
	- Update text and mark changes

- After supplements are added, double check all references in red to make sure they point to the right place