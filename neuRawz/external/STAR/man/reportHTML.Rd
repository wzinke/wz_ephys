\name{reportHTML}
\alias{reportHTML}
\title{Generic Function for Automatic HTML Report Generation}
\description{
  When a standard analysis is applied to some object it is useful to
  keep all the plots and summaries related to that analysis in a single
  place where they can be easily accessed and visualized. An html file
  containing the report of this analysis is ideally suited for that. The
  methods \code{reportHTML} generate such reports.
}
\usage{
reportHTML(object, filename, extension, directory, Title, ...)
}
\arguments{
  \item{object}{an object from which the report is going to be
    generated, perhaps following some standard analysis procedure.}
  \item{filename}{a character string. The generic name of all the files
    (html, png as well as \code{R} data files which will be
    generated. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{extension}{see \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{directory}{the full or relative path to the directory where the
    results are going to be stored. See also \code{\link[R2HTML]{HTMLInitFile}}.}
  \item{Title}{See \code{\link[R2HTML]{HTMLInitFile}}. If missing a default
    value baed on \code{filename} is provided.}
  \item{\dots}{additional parameters passed to the functions internally
    called by the actual methods.}
}
\value{
  Nothing is returned, an html file and figures in png format are
  written to disk together with the \code{R} variables generated during
  the analysis , if an analysis was performed.
}
\author{Christophe Pouzat \email{christophe.pouzat@gmail.com} }
\seealso{
  \code{\link{reportHTML.spikeTrain}},
  \code{\link{reportHTML.repeatedTrain}},
  \code{\link{reportHTML.gam}}
}
\examples{
##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
\keyword{file}
\keyword{IO}
