% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CytoProcessingStepImplementations.R
\name{applyFlowJoGate}
\alias{applyFlowJoGate}
\title{perform manual gating from a FlowJo gate}
\usage{
applyFlowJoGate(
  ff,
  gateName,
  wspSource = c("import", "pData"),
  wspFile = NULL,
  pDataVar = NULL,
  pDataPathMapping = NULL,
  ...
)
}
\arguments{
\item{ff}{a flowCore::flowFrame}

\item{gateName}{the name of the flowJo gate that will be applied}

\item{wspSource}{if "import", uses \code{wspFile} to read the flowjo workspace
(should be a .wsp file)
if "pData", uses \code{pDataVar} and \code{pDataPathMapping} to link a specific
phenotype data variable to map different flowjo workspace files}

\item{wspFile}{a flowjo workspace}

\item{pDataVar}{variable name (column of pheno data)
used to map the flowjo workspace file}

\item{pDataPathMapping}{a named list:
\itemize{
\item item names are possible values of \code{pDataVar}
\item item values are character() providing the \code{wspFile}
for the corresponding \code{pDataVar} value
}}

\item{...}{Extra arguments passed to \code{getFlowJoLabels()}}
}
\value{
a list with :
\itemize{
\item the first element ("matrix") is a matrix containing filtering results for
each specified gate
\item the second element ("labels") is a vector which assigns one label to each
cell. If no cell type correspond to a specific cell, than the keyword
'unlabeled' is assigned to this cell. If the cell belongs to several gates
(meaning that the gates ar not disjoints), than this cell is assigned to the
gate with the less matching cells.
}
}
\description{
perform manual gating by: reading a flowjo workspace file
using the flowWorkspace library, and apply one of the gates
to the inputflowFrame. Note that not only the selected gate, but also
all its hierarchy of gates, will be applied to the input dataset.
}
