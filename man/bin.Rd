% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin.R
\name{bin}
\alias{bin}
\title{Bin together data}
\usage{
bin(data, span, fx, .by = NULL)
}
\arguments{
\item{data}{Dataframe that contains the columns `time`, `id`, `x`, and `y`.}

\item{span}{Numeric denoting the size of the bins. Will pertain to the `time` 
variable, and should be specified in seconds (e.g., 0.5 means 500msec).}

\item{fx}{Function to execute on the data that falls within the bin. Should 
take in a similar data structure to `data` and output a single-row dataframe
with columns `x` and `y`.}

\item{.by}{String denoting whether the moving window should be taken with 
respect to a given grouping variable. Defaults to `NULL`}
}
\value{
Binned dataframe with a similar structure as `data`
}
\description{
Bin data together within a given window. While related, it is different from
`moving_window` in that the data quantity is reduced. Specifically, while 
`moving_window` keeps the same amount of data, `bin` will reduce the amount 
of data by taking all data within a bin together and summarizing it to one 
value.
}
\details{
A lot of this function is the same as for `moving window`: It only differs in 
its discrete jump to another bin instead of going from data point to data 
point.
}
