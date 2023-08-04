\name{sample_tau}
\alias{sample_tau}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample Kendall's tau
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample_tau(RVM, udata, STAN = STAN_model, iter_w = 150, iter_s = 2000, init = 0.8, s_bound = 0.95, refresh_s = 100, adapt_delta_s = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RVM}{
  A R-Vine Matrix object from the package \code{\link{VineCopula}}.
}
  \item{udata}{
  Copula data.
}
  \item{STAN}{
  A Stan model. (Default is a Stan model that estimates Kendall's tau for any R-Vine.)
}
  \item{iter_w}{
  Number of warmup iterations per Markov chain. (Default is \code{150}.)
}
  \item{iter_s}{
  Number of sampling iterations per Markov chain. (Default is \code{2000}.)
}
  \item{init}{
  Range or values for initial values for each Markov chain. (Default is \code{0.8}.)
}
  \item{s_bound}{
  Bound for Kendall's tau for all pair copulas in the R-Vine. If a pair copula can model positive and negative dependence, Kendall'S tau can be sampled in \code{[-s_bound, s_bound]}. If a pair copula can only model positive dependence, Kendall'S tau can be sampled in \code{[0, s_bound]}. If a pair copula can only model negative dependence, Kendall'S tau can be sampled in \code{[-s_bound, 0]}. (Default is \code{0.95}.)
}
  \item{refresh_s}{
  Number of iterations at which a message is printed. (Default is \code{100}.)
}
  \item{adapt_delta_s}{
  Target average acceptance probability. See Stan for more details. (Default is \code{0.8}.)
}
}
\details{
Jointly samples Kendall's tau for every edge in the given R-Vine. Pair copula families have to be known beforehand. For sampling, 4 Markov chains are used. Needs a valid \code{RVM} object as input as well as copula data. The \code{RVM} will be normalized if it isn't already.
}
\value{
A \code{\link{cmdstanr}} fit object containing the samples for all Kendall's tau and the corresponding pair copula parameters.
}
\references{
\url{https://mc-stan.org/docs/stan-users-guide/index.html} \cr
\url{https://cran.r-project.org/web/packages/VineCopula/VineCopula.pdf} \cr
\code{\link{cmdstanr}}
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
