setClass("blm",
               representation(
			       fit = "list",
			       par.start = "list",
			       f.loglik = "function",
			       f.score = "function",
			       run.time = "numeric",
			       data = "data.frame",
			       formula = "formula",
			       constraints = "list",
                               active.constraints = "list",
			       ineq = "matrix",
                               n.missing = "numeric",
                               H = "matrix",
                               V = "matrix"
			)
	)


setClass("lexpit",
               representation(
			       fit = "list",
			       par.start = "list",
			       f.loglik = "function",
			       f.score = "function",
			       run.time = "numeric",
			       data = "data.frame",
			       formula.linear = "formula",
                               formula.expit = "formula",
			       constraints = "list",
                               active.constraints = "list",
			       ineq = "list",
                               n.missing = "numeric",
                               H = "matrix",
                               V = "matrix"
			)
	)




