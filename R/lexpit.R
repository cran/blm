lexpit <- function(f.linear,f.expit,data,par.init,ineq=NULL,trace=FALSE,tol=1e-6,augmented=TRUE,warn=-1,...){

       warn.setting <- getOption("warn")
       options(warn=warn)
       
       if(missing(data)){
         stop("Dataset not supplied. Data must be given as a dataframe.")
       }

       if(missing(f.linear)|class(f.linear)!="formula"){
         stop("Invalid formula for linear component of model supplied.")
       }

       if(missing(f.expit)|class(f.expit)!="formula"){
         stop("Invalid formula for expit component of model supplied.")
       }
       
       na.rm <- remove.missing.lexpit(f.linear,f.expit,data)
       data = na.rm$data
       
       if(!is.na(na.rm$missing)){
         warning(paste(na.rm$missing,"rows with missing observation removed."))
       }
       else{
         na.rm$missing <- 0
       }
        
        if(missing(par.init)){
            par.start <- starting.values.lexpit(f.linear,f.expit,data)
        }
        else{
            par.start <- starting.values.lexpit(f.linear,f.expit,data,par.init=par.init)
         }


   LL <- lexpit.loglik(f.linear,f.expit,data)
   score <- lexpit.dot.loglik(f.linear,f.expit,data)

   constraints <- lexpit.constraints(f.linear,f.expit,data,ineq.mat=ineq)
       
   process.start = proc.time()

      if(augmented){
        
        fit <- auglag(par=par.start$par.start,fn=LL,gr=score,										hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
          }
        else{

       fit <- constrOptim.nl(par=par.start$par.start,fn=LL,gr=score,									hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
        }

   run.time = proc.time()-process.start
   run.time = as.numeric(run.time)[3]
   fit$par <- matrix(fit$par,ncol=1)
   rownames(fit$par) <- par.start$names
        
        	lexpit.object <- new("lexpit",
					  fit = fit,
					  par.start = par.start,
					  f.loglik = LL,
					  f.score = score,
					  run.time = run.time,
					  data = data,
					  formula.linear = f.linear,
					  formula.expit = f.expit,
					  constraints = constraints,
                                          active.constraints = list(),
       					  ineq = constraints$ineq.mat,
                                          n.missing=na.rm$missing,
                                          H = matrix(),
                                          V = matrix()
					  )

       if(augmented){
          H = -fit$hessian
          lexpit.object@active.constraints = check.auglag.lexpit.active.constraints(lexpit.object)
        }
       else{
          H = hessian.lexpit(lexpit.object)
          lexpit.object@active.constraints = check.lexpit.active.constraints(lexpit.object,tol)
        }
       
          V = -solve(H)

          lexpit.object@H = H
          lexpit.object@V = V
                
        options(warn=warn.setting)

        if(!augmented&!is.null(lexpit.object@active.constraints)){
          warning("\nEstimates at the boundary and augmented Lagrangian not used.\nStandard errors might be inaccurate.",immediate.=TRUE,call.=FALSE)
         }
       

 lexpit.object		
}
