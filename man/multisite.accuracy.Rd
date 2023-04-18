\name{multisite.accuracy}
\alias{multisite.accuracy}
\title{Estimate accuracy in multisite machine learning studies}
\description{Function to estimate the accuracy of the predictions of a multisite machine-learning model, controlling the effects of the site.}
\usage{
multisite.accuracy(y, y.pred, site, estimate = c("auc", "bac", "cor", "hr", "mse"),
                   site.method = "covar", mixed = FALSE, min.n = 10, ...)
}
\arguments{
  \item{y}{actual value of the variable that you want to predict. For \code{"auc"} estimates, it must be a binary vector. For \code{"bac"} estimates, it must be either a binary vector or a factor. For \code{"cor"} and \code{"mse"} estimates, it must be a numeric vector. For \code{"hr"} estimates, it must be an object of class \code{"Surv"}.}
  \item{y.pred}{predicted value of the variable. For \code{"auc"}, \code{"cor"}, \code{"hr"}, and \code{"mse"} estimates, it must be a numeric vector. For \code{"bac"} estimates, it must be either a binary vector or a factor with the same levels as \code{"y"}.}
  \item{site}{vector with the name of the site for each observation. Ignored if \code{site.method} is \code{"none"}.}
  \item{estimate}{\code{"auc"} (area under the ROC curve, when \code{"y"} is a binary variable and \code{"y.pred"} is a continous variable), \code{"bac"} (balanced accuracy, along with sensitivity and specificity when both \code{"y"} and \code{"y.pred"} are binary variables; or along with the sensitivity to detect each group when both \code{"y"} and \code{"y.pred"} are factors), \code{"cor"} and \code{"mse"} (correlation and mean squared error, when both \code{"y"} and \code{"y.pred"} are numeric variables), \code{"hr"} (hazard ratio, when \code{"y"} is an object of class \code{"Surv"} and \code{"y.pred"} is a numeric variable).}
  \item{site.method}{\code{"covar"} (site as covariate, e.g., using \code{lm}), \code{"meta"} (meta-analysis across sites), or \code{"none"} (no control of the effects of the site).}
  \item{mixed}{try to fit mixed-effects models, only for \code{site.method} = \code{"covar"} and \code{estimate} = \code{"cor"}, \code{"mse"}, or \code{"hr"}.}
  \item{min.n}{sites below this minimum sample size will be ignored.}
  \item{\dots}{further arguments for \code{logistf}, \code{lm}, \code{rma}, etc.}
}
\details{We refer the reader to the publication below for details on the calculations.}
\value{A data frame with the estimated accuracy, the methods used, and any warning or error.}
\references{
  Solanes, A., Palau, P., Fortea, L., Salvador, R., Gonzalez-Navarro, L., Llach, C.D., Valenti, M., Vieta, E., Radua, J. (2021)
  Biased accuracy in multisite machine-learning studies due to incomplete removal of the effects of the site.
  \emph{Psychiatry Research: Neuroimaging}, 314:111313.
  Solanes, A., Gosling, C.J., Fortea, L., Ortuno, M., Lopez-Soley, E., Llufriu, S., Madero, S., Martinez-Heras, E., Pomarol-Clotet, E., Solana, E., Vieta, E., Radua, J. (2021)
  Removing the effects of the site in brain imaging machine-learning - Measurement and extendable benchmark.
  \emph{Neuroimage}, 265:119800.
}
\examples{
for (real_effects in c(FALSE, TRUE)) {
  for (eos in c(FALSE, TRUE)) {
    
    # Simulate data
    site = rep(c("Site A", "Site B"), 200)
    y = c(scale(rnorm(400)))
    y.pred = c(scale(residuals(lm(rnorm(400) ~ y))))
    # If real effects:
    if (real_effects) {
      y.pred = c(scale(y.pred + y))
    }
    # If effects of the site:
    if (eos) {
      y = c(scale(y + 2 * (site == "Site B")))
      y.pred = c(scale(y.pred + 2 * (site == "Site B")))
    }
    cat("\n=== Real effects:", ifelse(real_effects, "yes", "no"),
        "\n=== Effects of the site:", ifelse(eos, "yes", "no"), "\n\n")
    
    # Numeric: without real effects, mse.pred should not be < mse.mean and cor should be ~0
    print(rbind(
      multisite.accuracy(y, y.pred, site, "mse", site.method = "covar"),
      multisite.accuracy(y, y.pred, site, "mse", site.method = "meta"),
      multisite.accuracy(y, y.pred, site, "mse", site.method = "none"))[,1:5])
    print(rbind(
      multisite.accuracy(y, y.pred, site, "cor", site.method = "covar"),
      multisite.accuracy(y, y.pred, site, "cor", site.method = "meta"),
      multisite.accuracy(y, y.pred, site, "cor", site.method = "none"))[,1:3])

  }
}
}
