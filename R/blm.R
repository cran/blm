blm <- function(f,data,par.init,ineq=NULL,trace=FALSE,tol=1e-6,adaptive=FALSE,...){

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
 
      LL <- blm.loglik(f,data)
      score <- blm.dot.loglik(f,data)
	

      constraints <- blm.constraints(f,data,A=ineq)
      
      process.start = proc.time()

      if(adaptive){
        
        fit <- auglag(par=par.start$par.start,fn=LL,gr=score,										hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
          }
        else{

       fit <- constrOptim.nl(par=par.start$par.start,fn=LL,gr=score,										hin=constraints$ineq,hin.jac=constraints$ineq.jac,
                                     control.outer=list(trace=trace,...))
        }

      run.time = proc.time()-process.start
      run.time = as.numeric(run.time)[3]
        
      if(is.null(ineq)) ineq = matrix()

      
		blm.object <- new("blm",
					  fit = fit,
					  par.start = par.start,
					  f.loglik = LL,
					  f.score = score,
					  run.time = run.time,
					  data = data,
					  formula = f,
					  constraints = constraints,
                                          active.constraints = list(),
					  ineq = constraints$A,
                                          n.missing = na.rm$missing,
                                          H = matrix(),
                                          V = matrix()
					  )

          H = hessian.blm(blm.object)
          V = -solve(H)

          blm.object@H = H
          blm.object@V = V
                
          blm.object@active.constraints = check.active.constraints(blm.object,tol)
        
        blm.object
}
