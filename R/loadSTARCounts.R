#' Load STAR aligner gene count output
#'
#' This function assumes that STAR aligner has been run with the
#' \code{--quantMode GeneCounts} option. It simply reads in the data
#' files, one-per-sample, into a matrix with rownames corresponding
#' to the gene names in the STAR output and column names that are
#' the file names passed to the function as the first parameter.
#'
#' The choice of \code{quant} is based on how STAR interprets the
#' reads relative to strandedness of the protocol used to produce
#' the data.  In particular, \code{unstranded} is used for an
#' unstranded protocol, \code{FRaligned} is for when the first read
#' is in the same direction as the strand of the mRNA, and
#' \code{SRaligned} is for when the second read is in the same
#' direction as the strand of the mRNA.
#'
#' @param fnames The file names for each of the .tab files from STAR.
#' @param quant Which column from the STAR output to use. See details section
#' for more information.
#' @return A matrix of counts, with genes from the STAR output on the rownames
#' and file names on the columns.
#'
#' @export
#' @importFrom readr read_tsv
#' @examples
#' fnames = dir(system.file(package='SeansStuff','extdata'),full.names=TRUE)
#' mat = loadSTARCounts(fnames)
#' dim(mat)
#' head(mat)
#' # clean up column names
#' colnames(mat) = basename(colnames(mat))
#' head(mat)
loadSTARCounts = function(fnames,quant=c('unstranded','FRaligned','SRaligned')) {
  choices = c('unstranded','FRaligned','SRaligned')
  col = match(match.arg(quant,choices),choices)+1
  geneNames = read_tsv(fnames[1],skip=4)[,1]
  mat = do.call('cbind',
                lapply(fnames,function(fname) {
                  read_tsv(fname,skip=4)[,col]
                  }))
  rownames(mat)=make.unique(geneNames)
  colnames(mat)=fnames
  return(mat)
}
