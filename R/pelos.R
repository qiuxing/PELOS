pelos <- function( COEFFICIENT , INITIAL.VALUE , INTERCEPT = NULL , OBSERVATION , tol = 1e-5 , max_iter = 30 , dif_step = 1e-3 , timestep = 1e-2 , Nprint = 10 )
#, sparse.output = FALSE )
{
    Nequation <- length(INITIAL.VALUE)
    Ntime <- length(TIMEPOINT)
    TIMEPOINT <- sort(unique(OBSERVATION[,2]))
#    NONZERO_ROW <- COEFFICIENT[,1]
#    NONZERO_COLUMN <- COEFFICIENT[,2]
#    COEFFICIENT.MATRIX <- COEFFICIENT[,3]
    Nnonzero <- length(COEFFICIENT)/3
    Nnonzero_intercept <- length(INTERCEPT)/2
    Nobservation <- length(OBSERVATION)/4
    max_func <- max_iter * ( Nequation+Nnonzero+Nnonzero_intercept+3 )
    Y <- numeric(Nobservation)

    res <- .Fortran ( 'S_lmdif_ODE_linear' , as.integer(Nequation) , as.integer(Ntime) , as.numeric(timestep) , as.numeric(TIMEPOINT) , as.integer(Nnonzero) , as.integer(COEFFICIENT[,1]) , as.integer(COEFFICIENT[,2]) , as.integer(Nnonzero_intercept) , as.integer(INTERCEPT[,1]) , as.integer(Nobservation) , as.integer(OBSERVATION[,1]) , as.numeric(OBSERVATION[,2]) , as.numeric(OBSERVATION[,3]) , as.numeric(OBSERVATION[,4]) , as.numeric(tol) , as.integer(max_func) , as.numeric(dif_step) , as.integer(Nprint) , as.numeric(c(INITIAL.VALUE,COEFFICIENT[,3],INTERCEPT[,2])) , as.numeric(Y) )
#    if ( sparse.output == FALSE )
#    {
        PM <- matrix(0,Nequation,Nequation)
        for ( i in 1:Nnonzero )
        {
            PM [ COEFFICIENT[i,1] , COEFFICIENT[i,2] ] <- res[[19]][Nequation+i]
        }
        IN <- numeric(Nequation)
        for ( i in 1:Nnonzero_intercept )
        {
            IN[INTERCEPT[i,1]] <- res[[19]][Nequation+Nnonzero+i]
        }
        return( list( initial.value = res[[19]][1:Nequation] , coefficient.matrix = PM , intercept = IN ) )
#    }
#    else
#    {
#        return( list( initial.value = res[[17]][1:Nequation] , coefficient.matrix = res[[17]][(Nequation+1):(Nequation+Nnonzero)] , intercept = res[[17]][(Nequation+Nnonzero+1):(Nequation+Nnonzero+Nnonzero_intercept)] ) )
#    }
}
