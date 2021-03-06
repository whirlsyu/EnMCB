% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univ_coxph.R
\name{univ_coxph}
\alias{univ_coxph}
\title{univariate and multivariate survival analysis using coxph}
\usage{
univ_coxph(dataframe, y_surv, digits = 4, asnumeric = TRUE)
}
\arguments{
\item{dataframe}{Clinic data and covariates ready to be tested. Rows are variables and columns are samples.}

\item{y_surv}{Survival function contain survival data, usually are obtained form Surv() function in survival package.}

\item{digits}{Integer indicating the number of decimal places.}

\item{asnumeric}{indicator that the data will be (True) / not (False) transformed into numeric. Default is true.}
}
\value{
Object of class \code{matrix} with results.
}
\description{
univariate and multivariate survival analysis using coxph
}
\examples{
data(demo_survival_data)
data('demo_data',package = "EnMCB")
demo_set<-demo_data$realdata
res<-univ_coxph(demo_set,demo_survival_data)
}
\author{
Xin Yu
}
