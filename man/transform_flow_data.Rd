% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{transform_flow_data}
\alias{transform_flow_data}
\title{Transform fluorescence values}
\usage{
transform_flow_data(flow_object, transform)
}
\arguments{
\item{flow_object}{A flow_object}

\item{transform}{Transformation method. Possible choices are 'asinh', 'none', 'log' or 'cytofAsinh'.}
}
\value{
A data.frame
}
\description{
Applies a transformation to fluorescence values embeeded in a Flow Object. Possible choices are hyperbolic sine (asinh)
using cofactors inferred by `find_cofactors()`, biexponential (none), logarithmic (in base 10), or CYTOF asinh (cofactor of 5).
}
