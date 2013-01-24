\name{pelos}
\alias{pelos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parameter Estimation of Linear ODE System }
\description{ An optimization algorithm based on Levenberg-Marquardt algorithm to minimize the total sum of squares between observation value and numerical solution of the ODE system.  It returns the optimum initial value, intercept and coefficient matrix according to the sparse structure given. }
\usage{
pelos (TIMEPOINT, COEFFICIENT, INITIAL.VALUE, INTERCEPT, OBSERVATION,
tol=1e-5, max_iter=30, dif_step=1e-3, timestep = 1e-2,
Nprint=10, sparse.output=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{TIMEPOINT}{A vector which gives time-points of observations.}
    \item{COEFFICIENT}{
        A 3-column matrix whcih gives sparse structure and starting value of coefficient matrix of linear ODE system.  Each row stands for one non-zero entry.  The first two columns are row,column number respectively, and the third column is its starting value.
    }
    \item{INITIAL.VALUE,INTERCEPT}{
        Starting value of initial value and intercept of the ODE system respectively.  Must be of same length.
    }
    \item{OBSERVATION}{
        A 4-column matrix which gives observation data used for sum-of-squares object fucntion.  Each row stands for one observation point.  The first three columns are its index,time-point and value respectively.  The last column is its weight in the object function.  Notice that this weight is inside the square.
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
    \item{sparse.output}{
        If true, the output coefficient matrix in return value is in sparse format, which is a vector of the same order as the third column of input parameter COEFFICIENT.  Otherwise, the return value is a full square matrix with 0 filled in.
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
%\examples{
%a <- 1:5
%b <- 5:20
%y <- conv(a,b)
%y
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{least-squares}
