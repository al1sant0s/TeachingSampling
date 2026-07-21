#' @title Full Business Population Database
#' @description
#' A data frame containing socioeconomic and financial variables from a
#' population of 2,396 business units, used throughout the package to
#' illustrate survey sampling designs and estimators.
#' @format A data frame with 2,396 rows and 8 variables:
#' \describe{
#'   \item{ID}{The identifier of the company. It corresponds to an
#'     alphanumeric sequence (two letters and ten digits).}
#'   \item{Ubication}{The address of the principal office of the company.}
#'   \item{Level}{The size level of the company discriminated according to
#'     the income declared. There are small, medium and big companies.}
#'   \item{Zone}{The country is divided by counties. A company belongs to a
#'     particular zone according to its cartographic location.}
#'   \item{Income}{The total amount of a company's earnings in the previous
#'     fiscal year.}
#'   \item{Employees}{The total number of persons working for the company in
#'     the previous fiscal year.}
#'   \item{Taxes}{The total amount of a company's income tax.}
#'   \item{SPAM}{Indicates if the company uses the Internet and webmail
#'     options to make self-propaganda.}
#' }
#' @seealso \code{\link{BigLucy}}, \code{\link{BigCity}}
#' @author Hugo Andres Gutierrez Rojas \email{hagutierrezro@@gmail.com}
#' @references
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' @examples
#' data(Lucy)
#' attach(Lucy)
#' # The variables of interest are: Income, Employees and Taxes
#' estima <- data.frame(Income, Employees, Taxes)
#' # The population totals
#' colSums(estima)
#' # Some parameters of interest
#' table(SPAM, Level)
#' xtabs(Income ~ Level + SPAM)
#' # Correlations among characteristics of interest
#' cor(estima)
#' # Some useful histograms
#' hist(Income)
#' hist(Taxes)
#' hist(Employees)
#' # Some useful plots
#' boxplot(Income ~ Level)
#' barplot(table(Level))
#' pie(table(SPAM))
"Lucy"

#' @title Full Business Population Database (Extended)
#' @description
#' A data frame corresponding to some financial variables of 85,396 industrial
#' companies of a city in a particular fiscal year.
#' @format A data frame with the following variables:
#' \describe{
#'   \item{ID}{The identifier of the company. It corresponds to an
#'     alphanumeric sequence (two letters and three digits).}
#'   \item{Ubication}{The address of the principal office of the company
#'     in the city.}
#'   \item{Level}{The industrial companies are discriminated according to the
#'     income declared. There are small, medium and big companies.}
#'   \item{Zone}{The country is divided by counties. A company belongs to a
#'     particular zone according to its cartographic location.}
#'   \item{Income}{The total amount of a company's earnings in the previous
#'     fiscal year.}
#'   \item{Employees}{The total number of persons working for the company in
#'     the previous fiscal year.}
#'   \item{Taxes}{The total amount of a company's income tax.}
#'   \item{SPAM}{Indicates if the company uses the Internet and webmail
#'     options to make self-propaganda.}
#'   \item{ISO}{Indicates if the company is certified by the International
#'     Organization for Standardization.}
#'   \item{Years}{The age of the company.}
#'   \item{Segments}{Cartographic segments by county. A segment comprises
#'     on average 10 companies located close to each other.}
#' }
#' @seealso \code{\link{Lucy}}, \code{\link{BigCity}}
#' @author Hugo Andres Gutierrez Rojas \email{hagutierrezro@@gmail.com}
#' @references
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' @examples
#' data(BigLucy)
#' attach(BigLucy)
#' # The variables of interest are: Income, Employees and Taxes
#' estima <- data.frame(Income, Employees, Taxes)
#' # The population totals
#' colSums(estima)
#' # Some parameters of interest
#' table(SPAM, Level)
#' xtabs(Income ~ Level + SPAM)
#' # Correlations among characteristics of interest
#' cor(estima)
#' # Some useful histograms
#' hist(Income)
#' hist(Taxes)
#' hist(Employees)
#' # Some useful plots
#' boxplot(Income ~ Level)
#' barplot(table(Level))
#' pie(table(SPAM))
"BigLucy"

#' @title Full Person-level Population Database
#' @description
#' A data frame corresponding to some socioeconomic variables from 150,266
#' people of a city in a particular year.
#' @format A data frame with the following variables:
#' \describe{
#'   \item{HHID}{The identifier of the household. It corresponds to an
#'     alphanumeric sequence (four letters and five digits).}
#'   \item{PersonID}{The identifier of the person within the household.
#'     Note: it is not a unique identifier for the whole population.}
#'   \item{Stratum}{Households are located in geographic strata. There are
#'     119 strata across the city.}
#'   \item{PSU}{Households are clustered in cartographic segments defined as
#'     primary sampling units (PSU). There are 1,664 PSU nested within strata.}
#'   \item{Zone}{Segments within strata can be located in urban or rural
#'     areas across the city.}
#'   \item{Sex}{Sex of the person.}
#'   \item{Age}{Age of the person.}
#'   \item{Income}{Per capita monthly income.}
#'   \item{Expenditure}{Per capita monthly expenditure.}
#'   \item{Employment}{A person's employment status.}
#'   \item{Poverty}{Indicates whether the person is poor or not, based on
#'     income.}
#'   \item{MaritalST}{Marital status of the person.}
#' }
#' @seealso \code{\link{Lucy}}, \code{\link{BigLucy}}
#' @author Hugo Andres Gutierrez Rojas \email{hagutierrezro@@gmail.com}
#' @references
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas
#' y estimacion de parametros}. Editorial Universidad Santo Tomas.
#' @examples
#' data(BigCity)
#' attach(BigCity)
#' estima <- data.frame(Income, Expenditure)
#' # The population totals
#' colSums(estima)
#' # Some parameters of interest
#' table(Poverty, Zone)
#' xtabs(Income ~ Poverty + Zone)
#' # Correlations among characteristics of interest
#' cor(estima)
#' # Some useful histograms
#' hist(Income)
#' hist(Expenditure)
#' # Some useful plots
#' boxplot(Income ~ Poverty)
#' barplot(table(Employment))
#' pie(table(MaritalST))
"BigCity"