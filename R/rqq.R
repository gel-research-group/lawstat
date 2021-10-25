#' Test of Normality Using RQQ Plots
#' 
#' Produce robust quantile-quantile (RQQ) and classical quantile-quantile (QQ) 
#' plots for graphical assessment of normality and optionally add a line, a QQ line, 
#' to the produced plot. The QQ line may be chosen to be a 45-degree line or to pass
#' through the first and third quartiles of the data. 
#' \code{NA}s from the data are omitted.
#' 
#' @details An RQQ plot is a modified QQ plot where data are robustly standardized 
#' by the median and robust measure of spread (rather than mean and classical 
#' standard deviation as in the basic QQ plots) and then are plotted against the 
#' expected standard normal order statistics 
#' \insertCite{Gel_etal_2005,Weisberg_2005}{lawstat}. 
#' Under normality, the plot of the standardized 
#' observations should follow the 45-degree line, or QQ line. Both the median and robust 
#' standard deviation are significantly less sensitive to outliers than mean and 
#' classical standard deviation and therefore are more preferable in many practical 
#' situations to assess graphically deviations from normality (if any). We choose 
#' median and MAD as a robust measure of location and spread for our RQQ plots since 
#' this standardization typically provides a clearer graphical diagnostics of normality. 
#' In particular, deviations from the QQ line are usually more noticeable in RQQ plots 
#' in the case of outliers and heavy tails. Users can also choose to plot the 
#' 45-degree line or the 1st-3rd quartile line (see the argument \code{line.type}). 
#' No line is the default.
#'
#'
#' @param y the input data.
#' @param plot.it logical. Should the result be plotted?
#' @param square.it logical. Should the plot scales be square? The default is \code{TRUE}.
#' @param scale the choice of a scale estimator, i.e., the classical or robust estimate 
#' of the standard deviation.
#' @param location the choice of a location estimator, i.e., the mean or median.
#' @param line.it logical. Should the line be plotted? No line is the default.
#' @param line.type If \code{line.it = TRUE}, the choice of a line to be plotted, i.e., 
#' the 45-degree line or the line passing through the first and third quartiles 
#' of the data.
#' @param col.line the color of the line (if plotted).
#' @param lwd the line width (if plotted).
#' @param outliers logical. Should the outliers be listed in the output?
#' @param alpha significance level of outliers. If \code{outliers = TRUE}, then all 
#' observations that are less than the \code{100*alpha}-th standard normal percentile or 
#' greater than the \code{100*(1-alpha)}-th standard normal percentile will be listed 
#' in the output.
#' @param ... other parameters passed to the \code{\link[graphics]{plot}} function.
#'
#'
#' @return A list with the following numeric components:
#' \item{x}{the x-coordinates of the points that were/would be plotted.}
#' \item{y}{the original data vector, i.e., the corresponding y-coordinates,
#' including \code{NA}s (if any).}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{rjb.test}}, \code{\link{sj.test}}, 
#' \code{\link[stats]{qqnorm}}, \code{\link[stats]{qqplot}}, \code{\link[stats]{qqline}}
#' 
#' @keywords distribution robust
#' 
#' @author W. Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' ## Simulate 100 observations from standard normal distribution:
#' y = rnorm(100)
#' rqq(y)
#' 
#' ## Using Michigan data
#' data(michigan)
#' rqq(michigan)
#' 
`rqq` <- function (y,
                   plot.it = TRUE,
                   square.it = TRUE,
                   scale = c("MAD", "J", "classical"),
                   location = c("median", "mean"),
                   line.it = FALSE,
                   line.type = c("45 degrees", "QQ"),
                   col.line = 1,
                   lwd = 1,
                   outliers = FALSE,
                   alpha = 0.05,
                   ...)
{
    y <- na.omit(y)
    x = sort(y)
    scale <- match.arg(scale)
    location <- match.arg(location)
    line.type <- match.arg(line.type)
    if (location == "mean") {
        M = mean(x)
    }
    else {
        M = median(x)
    }
    if (scale == "classical") {
        qqstd = "QQ plot standardized by the classical std dev and"
        y = (x - M) / sd(x)
    }
    else if (scale == "MAD") {
        qqstd = "RQQ plot standardized by MAD and"
        y = (x - M) / mad(x)
    }
    else {
        scale = "J"
        qqstd = "RQQ plot standardized by J and"
        j = sqrt(pi / 2) * mean(abs(x - median(x)))
        y = (x - M) / j
    }
    if ((line.it == "TRUE") & (line.type == "QQ")) {
        qql = ", QQ line"
    }
    else if ((line.it == "TRUE") & (line.type == "45 degrees")) {
        qql = ", 45 degrees line"
    }
    
    qq <- qqnorm(y, plot.it = FALSE)
    if (line.it == TRUE) {
        if (square.it == TRUE) {
            q <- qqnorm(
                y,
                xlim = c(min(qq$x, y) - 0.5, max(qq$x, y) + 0.5),
                ylim = c(min(qq$x, y) - 0.5, max(qq$x, y) + 0.5),
                main = paste(qqstd, location, qql),
                ...
            )
        }
        else {
            q <- qqnorm(y, main = paste(qqstd, location, qql),
                        ...)
        }
        if (line.type == "QQ") {
            qqline(y,
                   datax = FALSE,
                   col = col.line,
                   lwd = lwd)
        }
        else {
            abline(0, 1, col = col.line, lwd = lwd)
        }
    }
    else {
        if (square.it == TRUE) {
            q <- qqnorm(
                y,
                xlim = c(min(qq$x, y) - 0.5, max(qq$x, y) + 0.5),
                ylim = c(min(qq$x, y) - 0.5, max(qq$x, y) + 0.5),
                main = paste(qqstd, location),
                ...
            )
            
        }
        else {
            q <- qqnorm(y, main = paste(qqstd, location), ...)
        }
    }
    
    if (outliers == TRUE) {
        left.out = y[y < qnorm(alpha)]
        right.out = y[y > qnorm(1 - alpha)]
        print(data.frame(left.tail.outliers = left.out))
        print(data.frame(right.tail.outliers = right.out))
    }
    
}
