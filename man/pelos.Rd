\name{pelos}
\alias{pelos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parameter Estimation of Linear ODE System }
\description{ An optimization algorithm based on Levenberg-Marquardt algorithm to minimize the total sum of squares between observation value and numerical solution of the ODE system.  It determines the initial value, intercept and coefficient matrix according to the sparse structure given. }
\usage{
pelos ( COEFFICIENT , INITIAL.VALUE , INTERCEPT = NULL , OBSERVATION , tol = 1e-5 , max_iter = 30 , dif_step = 1e-3 , timestep = 1e-2 , Nprint = 10 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{COEFFICIENT}{
        A 3-column matrix which gives sparse structure and starting value of coefficient matrix of linear ODE system.  Each row stands for one non-zero entry.  The first two columns are row,column number respectively, and the third column is its starting value.
    }
    \item{INITIAL.VALUE}{
        Starting value of initial value.  Its length must be equal to number of equations.
    }
    \item{INTERCEPT}{
        A 2-column matrix which gives sparse structure and starting value of intercept vector.  Each row stands for one non-zero entry.  The first column is position of the non-zero entry and second column is its starting value.  Its default value is NULL, which forces intercept to be a 0-vector.
    }
    \item{OBSERVATION}{
        A 4-column matrix which gives observation data used for sum-of-squares object function.  Each row stands for one observation point.  The first three columns are its index,time-point and value respectively.  The last column is its weight in the object function.  Notice that this weight is inside the square.
    }
    \item{tol}{
        Tolerance parameter to terminate iteration.  A smaller value gives a more accurate solution but allows the algorithm to perform more iterations.
    }
    \item{max_iter}{
        Algorithm is forced to terminate after this many iterations.
    }
    \item{dif_step}{
        Step length to perform finite difference.  This value should be larger than noise level.
    }
    \item{timestep}{
        Time step length to solve ODE system in Runge-Kutta algorithm.  A smaller value gives a more accurate ODE numerical solution but takes more time in each iteration.
    }
    \item{Nprint}{
        Output sum of squares value after every this many iterations.  Setting it to 0 to disable outputs in iterative process.
    }
}
%\details{
%  More details about this function.  May be very long (but please don't
%  do that!).
%}
\value{
    A list of three items: initial.value,intercept,coefficient.matrix.
}
%\references{
%  Blah blah blah.
%}
\author{
    GuuD. <guuuuud@gmail.com>
}
%\note{
%  blah blah blah
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%  \code{\link{sum}}, \code{\link{help}}.
%}
\examples{
data(Cells)
res <- pelos(A0, X0, b0, Cells, max_iter=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{least-squares}
