\name{SysCall}
\alias{SysCall}
\title{Wrapper for system commands}
\usage{
SysCall(exec, ..., args = list(), stdin = NULL, stdout = NULL,
  redirection = TRUE, style = c("unix", "gnu"), sep = " ",
  show_cmd = FALSE, intern = FALSE, input = NULL)
}
\arguments{
  \item{exec}{The system command to be invoked.}

  \item{...}{Arguments passed on to the \code{system}
  command as name-value or name=\code{TRUE} pairs.}

  \item{args}{Named list of arguments passed on to the
  \code{system} command. Is merged with \code{...}.}

  \item{stdin}{Input.}

  \item{stdout}{Output.}

  \item{redirection}{Redirection.}

  \item{style}{One of \sQuote{unix} or \sQuote{gnu}.}

  \item{sep}{Seperator of option and option argument.}

  \item{show_cmd}{Have a look what the final command looks
  like.}

  \item{intern}{Passed on to \code{\link{system}}'s
  \code{intern} argument.}

  \item{input}{Passed on to \code{\link{system}}'s
  \code{input} argument.}
}
\description{
Wrapper for system commands
}
\keyword{internal}

