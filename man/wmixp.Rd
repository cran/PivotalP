\name{wmixp}
\alias{wmixp}
\alias{print.wmixp}
\title{Prediction future points from mixture weibull distribution}

\description{
Construct a prediction point for future observations from mixture weibull distribution. Generic method is \code{print}.
}
\usage{
wmixp(data, s, n, a ,parameters, conf=0.95)

\method{print}{wmixp}(x, \dots)
}

\arguments{
\item{data}{A numeric vector.}
\item{s}{
A numeric value the order of prediction point.}
\item{n}{
A numeric vector for the size of all data.}
\item{a}{
A numeric value of mixing proportion.}
\item{parameters}{ A numeric vector of the parameter of distributions}
\item{conf}{Confidence level for the test.}
\item{x}{An object of class \code{"wmixp"}.  }
\item{\dots}{Further argument to be passed to generic function}
}
\details{
Prediction of future observations if the data follows a mixture of two weibull distributions}
\value{
\code{wmixp} returns an object of class \code{"wmixp"}, a list with the following components:
\item{ interval }{ the prediction interval. }
\item{ lower }{ the lower bound of the interval. }
\item{ upper }{ the upper bound of the interval. }
\item{ r }{ the length of the data. }
\item{ s }{ the order of the next observation. }
\item{ n }{ the length of all the data. }
\item{ parameters }{ the parameter estimate. }

Generic function:
\describe{
\item{\code{print}}{
    The print of a \code{"wmixp"} object shows the prediction point(s) for the future observation(s).
}}
}

\seealso{ \code{PredictionR}.}

\references{

H. M. Barakat, Magdy E. El-Adll, Amany E. Aly (2014), \emph{Prediction intervals of future observations for a sample random size from any continuous distribution}.
Mathematics and Computers in Simulation, volume 97, 1-13.

O. M. Khaled, K. S. Khalil and M. H. Harby (2023), \emph{PREDICTING FUTURE DATA FROM GAMMA-MIXTURE AND BETA-MIXTURE DISTRIBUTIONS AND APPLICATION TO THE RECOVERY RATE OF COVID-19}.
Advances and Applications in Statistics (AAIS), OCT, 2023.
}
\author{
O. M. Khaled, K. S. Khalil and M. H. Harby.
}

\examples{
# prediction interval and point for the next observations based on mixture weibull distribution
#
set.seed(123)
x1 <- 0.5*rweibull(7, 4, 2)+0.5*rweibull(7, 1, 3)
wmixp(x1, 8, 10,0.5,c(4,2,1,3),conf=0.95)
}

