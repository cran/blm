print.like.glm <- function(x,digits=4){

  ## MODIFIED FROM print.glm
     beta = coef(x)
     p = nrow(beta)
     df = nrow(x@data)-p
     
     if(class(x)=="blm"){
     cat("\nCall:  ", paste(deparse(x@formula), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
   }
     else{
     cat("\nLinear Call:  ", paste(deparse(x@formula.linear), sep = "\n", collapse = "\n"),
         sep = "")
          cat("\nExpit Call:  ", paste(deparse(x@formula.expit), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
     }
        cat("Coefficients")
        cat(":\n")

        T <- beta/sqrt(diag(vcov(x)))
        V <- sqrt(diag(vcov(x)))
        P <- mapply(pt,q=abs(T),MoreArgs=list(df=df))
        P <- 2*(1-P)
        coeftable <- cbind(beta,T,V,P)
        coeftable <- format(coeftable,digits=digits,sci=TRUE)
        colnames(coeftable) <- c("estimates","t-value","std. err","p-value")
                       
        print.default(coeftable, 
        print.gap = 2, quote = FALSE)

    cat("\nDegrees of Freedom:",df,"\n")
    cat("Run time (sec):",x@run.time,"\n")
    LL <- -x@f.loglik(beta)
    AIC <- 2*length(beta)-2*LL
    cat("LogLik:\t   ", format(signif(LL, 
        digits)), "\tAIC:", format(signif(AIC, digits)))
                 
    invisible(x)
}


setMethod("print","blm",print.like.glm)
setMethod("print","lexpit",print.like.glm)
