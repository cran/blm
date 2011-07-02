ssc <- function(object){

    J <- qr(object@active.constraints$active.grad)
    N <- Null(object@active.constraints$active.grad)
    positive.definite <- all(eigen(-t(N)%*%object@H%*%N)$values>0)
    barrier <- object@fit$barrier.value

    list(barrier.value=barrier,is.positive.definite=positive.definite,row.rank=J$rank)
}
