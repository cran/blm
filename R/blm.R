blm <- function(f,data,par.init,weights=NULL,ineq=NULL,trace=FALSE,tol=1e-6,
                augmented=TRUE,warn=-1,...){

       warn.setting <- getOption("warn")
       options(warn=warn)
         
       if(missing(data)){
         stop("Dataset not supplied. Data must be given as a dataframe.")
       }

       if(missing(f)|class(f)!="formula"){
         stop("Invalid formula supplied.")
       }
  
       na.rm <- remove.missing.blm(f,data)
       data = na.rm$data

       if(!is.na(na.rm$missing)){
         warning(paste(na.rm$missing,"rows with missing observation removed."))
       }
       else{
         na.rm$missing <- 0
       }
       
       
        if(missing(par.init)){
           par.start <- starting.values.blm(f,data)
        }
        else{
           par.start <- starting.values.blm(f,data,par.init=par.init)
        }

      weighted = ifelse(is.null(weights),FALSE,TRUE)

      if(!weighted) weights = rep(1,nrow(data))
       
      LL <- blm.loglik(f,data,weights)
      score <- blm.dot.loglik(f,data,weights)
      constraints <- blm.constraints(f,data,ineq.mat=ineq)

      process.start = proc.time()

      if(augmented){
        
        fit <- auglag(par=par.start$par.start,fn=LL,gr=score,										hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
          }
        else{

       fit <- constrOptim.nl(par=par.start$par.start,fn=LL,gr=score,									 hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
        }

      run.time = proc.time()-process.start
      run.time = as.numeric(run.time)[3]
      fit$par <- matrix(fit$par,ncol=1)
      rownames(fit$par) <- par.start$names

		blm.object <- new("blm",
					  fit = fit,
					  par.start = par.start,
					  f.loglik = LL,
					  f.score = score,
                                          weights = weights,
                                          run.time = run.time,
					  data = data,
					  formula = f,
					  constraints = constraints,
                                          active.constraints = list(),
					  ineq = constraints$ineq.mat,
                                          n.missing = na.rm$missing,
                                          H = matrix(),
                                          V = matrix()
					  )

             
       if(augmented){
          H = -fit$hessian
          blm.object@active.constraints = check.auglag.blm.active.constraints(blm.object)
        }
       else{
          H = hessian.blm(blm.object)
          blm.object@active.constraints = check.active.constraints(blm.object,tol)
        }

          V = -solve(H)

          blm.object@H = H
          blm.object@V = V

       if(weighted){    #SANDWHICH ESTIMATOR
         S <- weighted.vcov.blm(blm.object)
         blm.object@V = V%*%S%*%V
       }
       
       options(warn=warn.setting)

        if(!augmented&!is.null(blm.object@active.constraints)){
          warning("\nEstimates at the boundary and augmented Lagrangian not used.\nStandard errors might be inaccurate.",immediate.=TRUE,call.=FALSE)
         }
       
 blm.object
}
