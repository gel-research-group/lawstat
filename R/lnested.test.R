#' Test for a Monotonic Trend in Variances
#'
#' The test statistic is based on the finite intersection approach.
#' 
#' @details The test statistic is based on 
#' the classical Levene's procedure (using the group means), 
#' the modified Brown--Forsythe Levene-type procedure (using the group medians),
#' or the modified Levene-type procedure (using the group trimmed means).
#' More robust versions of the test using the correction factor or structural zero
#' removal method are also available. Two options for calculating critical values,
#' namely, approximated and bootstrapped, are available.
#' By default, \code{NA}s are omitted from the data.
#'
#'
#' @param y a numeric vector of data values.
#' @param group factor of the data.
#' @param location the default option is \code{"median"} corresponding to the
#' robust Brown--Forsythe Levene-type procedure \insertCite{Brown_Forsythe_1974}{lawstat}; 
#' \code{"mean"} corresponds to the classical Levene's procedure
#' \insertCite{Levene_1960}{lawstat}, and \code{"trim.mean"} corresponds to the robust
#' Levene-type procedure using the group trimmed means.
#' @param tail the default option is \code{"right"}, corresponding to an increasing
#' trend in variances as the one-sided alternative;  \code{"left"} corresponds to a
#' decreasing trend in variances, and  \code{"both"} corresponds to any
#' (increasing or decreasing) monotonic trend in variances as the two-sided alternative.
#' @param trim.alpha the fraction (0 to 0.5) of observations to be trimmed from
#' each end of \code{x} before the mean is computed.
#' @param bootstrap a logical value identifying whether to implement bootstrap.
#' The default is \code{FALSE}, i.e., no bootstrap; if set to \code{TRUE},
#' the bootstrap method described in \insertCite{Lim_Loh_1996;textual}{lawstat}
#' for Levene's test is applied.
#' @param num.bootstrap number of bootstrap samples to be drawn when the \code{bootstrap}
#' argument is set to \code{TRUE}. The default value is 1000.
#' @param correction.method procedures to make the test more robust;
#' the default option is \code{"none"}; \code{"correction.factor"} applies the
#' correction factor described by \insertCite{OBrien_1978;textual}{lawstat} and
#' \insertCite{Keyes_Levy_1997;textual}{lawstat}; \code{"zero.removal"} performs the
#' structural zero removal method by \insertCite{Hines_Hines_2000;textual}{lawstat};
#' \code{"zero.correction"} performs a combination of the O'Brien's correction factor
#' and the Hines--Hines structural zero removal method \insertCite{Noguchi_Gel_2010}{lawstat}.
#' Note that the options \code{"zero.removal"} and \code{"zero.correction"} are only
#' applicable when the location is set to \code{"median"}, otherwise, \code{"none"} is applied.
#' @param correlation.method measures of correlation; the default option is
#' \code{"pearson"}, the linear correlation coefficient that is equivalent to the t-test;
#' nonparametric measures of correlation such as \code{"kendall"} (Kendall's tau)
#' or \code{"spearman"} (Spearman's rho) may also be chosen.
#'
#'
#' @return A list with the following elements:
#' \item{T}{the statistic and \eqn{p}-value of the test based on the Tippett \eqn{p}-value combination.}
#' \item{F}{the statistic and \eqn{p}-value of the test based on the Fisher \eqn{p}-value combination.}
#' \item{N}{the statistic and \eqn{p}-value of the test based on the Liptak \eqn{p}-value combination.}
#' \item{L}{the statistic and \eqn{p}-value of the test based on the Mudholkar--George \eqn{p}-value combination.}
#'
#' Each of the list elements is a list of class \code{"htest"} with the following elements:
#' \item{statistic}{the value of the test statistic expressed in terms of correlation (Pearson, Kendall, or Spearman).}
#' \item{p.value}{the \eqn{p}-value of the test.}
#' \item{method}{type of test performed.}
#' \item{data.name}{a character string giving the name of the data.}
#' \item{non.bootstrap.statistic}{the statistic of the test without bootstrap method.}
#' \item{non.bootstrap.p.value}{the \eqn{p}-value of the test without bootstrap method.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{levene.test}}, \code{\link{ltrend.test}},
#' \code{\link{mma.test}}, \code{\link{neuhauser.hothorn.test}},
#' \code{\link{robust.mmm.test}}
#'
#' @keywords htest robust variability
#'
#' @author Kimihiro Noguchi, W. Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#'
#' @export
#' @examples
#' data(pot)
#' lnested.test(pot[,"obs"], pot[, "type"], location = "median", tail = "left",
#'              correction.method = "zero.correction")$N
#'
#' lnested.test(pot[, "obs"], pot[, "type"], location = "median", tail = "left",
#'              correction.method = "zero.correction",
#'              bootstrap = TRUE, num.bootstrap = 500)$N
#'
`lnested.test` <- function(y,
                           group,
                           location = c("median", "mean", "trim.mean"),
                           tail = c("right", "left", "both"),
                           trim.alpha = 0.25,
                           bootstrap = FALSE,
                           num.bootstrap = 1000,
                           correction.method = c("none",
                                                 "correction.factor",
                                                 "zero.removal",
                                                 "zero.correction"),
                           correlation.method = c("pearson", "kendall", "spearman"))
{
    ### assign location, tail, and a correction method ###
    
    location <- match.arg(location)
    tail <- match.arg(tail)
    correlation.method <- match.arg(correlation.method)
    correction.method <- match.arg(correction.method)
    DNAME <- deparse(substitute(y))
    y <- y[!is.na(y)]
    group <- group[!is.na(y)]
    
    ### stop the code if the length of y does not match the length of group ###
    
    if (length(y) != length(group))
    {
        stop("the length of the data (y) does not match the length of the group")
    }
    
    ### sort the order just in case the input is not sorted by group ###
    
    reorder <- order(group)
    group <- group[reorder]
    y <- y[reorder]
    gr <- group
    group <- as.factor(group)
    
    ### define the measure of central tendency (mean, median, trimmed mean) ###
    
    if (location == "mean") {
        means <- tapply(y, group, mean)
        METHOD <-
            "lnested test based on classical Levene's procedure using the group means"
    } else if (location == "median") {
        means <- tapply(y, group, median)
        METHOD <-
            "lnested test based on the modified Brown-Forsythe Levene-type procedure using the group medians"
    } else {
        location = "trim.mean"
        means <- tapply(y, group, mean, trim = trim.alpha)
        METHOD <-
            "ltrend test based on the modified Brown-Forsythe Levene-type procedure using the group trimmed means"
    }
    
    ### calculate the sample size of each group and deviation from center ###
    z <- y - means[group]
    n <- tapply(z, group, length)
    ngroup <- n[group]
    
    ### multiply the correction factor to each observation if "correction.factor" is chosen ###
    
    if (correction.method == "correction.factor")
    {
        METHOD <- paste(METHOD, "with correction factor applied")
        correction <- sqrt(ngroup / (ngroup - 1))
        z <- z * correction
    }
    
    resp.mean <- abs(z)
    
    k <- length(n)
    m <- k - 1
    
    ### create vectors for storing component p-values (bootstrap and non-bootstrap) ###
    
    p <- double(m)
    q <- double(m)
    non.bootstrap.p <- double(m)
    non.bootstrap.q <- double(m)
    
    ### assign no correction technique if the central tendency is median, and ###
    ### any technique other than "correction.factor" is chosen                ###
    
    if (location != "median" &&
        correction.method != "correction.factor")
    {
        METHOD <-
            paste(
                METHOD,
                "(",
                correction.method,
                "not applied because the location is not set to median",
                ")"
            )
        correction.method = "none"
    }
    
    ### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
    ### or "zero.correction" (Noguchi and Gel, 2009).                            ###
    
    if (correction.method == "zero.correction" ||
        correction.method == "zero.removal")
    {
        if (correction.method == "zero.removal")
        {
            METHOD <-
                paste(METHOD,
                      "with Hines-Hines structural zero removal method")
        }
        if (correction.method == "zero.correction")
        {
            METHOD <-
                paste(
                    METHOD,
                    "with modified structural zero removal method and correction factor"
                )
        }
        
        ### set up variables for calculating the deviation from center ###
        
        resp.mean <- z
        k <- length(n)
        temp <- double()
        endpos <- double()
        startpos <- double()
        
        ### calculate the absolute deviation from mean and remove zeros ###
        
        for (i in 1:k)
        {
            group.size <- n[i]
            j <- i - 1
            
            ### calculate the starting and ending index of each group ###
            
            if (i == 1)
                start <- 1
            else
                start <- sum(n[1:j]) + 1
            startpos <- c(startpos, start)
            end <- sum(n[1:i])
            endpos <- c(endpos, end)
            
            ### extract the deviation from center for the ith group ###
            
            sub.resp.mean <- resp.mean[start:end]
            sub.resp.mean <- sub.resp.mean[order(sub.resp.mean)]
            
            ### remove structural zero for the odd-sized group ###
            
            if (group.size %% 2 == 1)
            {
                mid <- (group.size + 1) / 2
                temp2 <- sub.resp.mean[-mid]
            }
            ### remove structural zero for the even-sized group ###
            
            if (group.size %% 2 == 0)
            {
                mid <- group.size / 2
                
                ### set up the denominator value for the transformation ###
                
                ### set 1 for the "zero.correction" option ###
                
                denom <- 1
                
                ### set sqrt(2) for the "zero.removal" option ###
                
                if (correction.method == "zero.removal")
                    denom <- sqrt(2)
                
                ### perform the orthogonal transformation ###
                
                replace1 <-
                    (sub.resp.mean[mid + 1] - sub.resp.mean[mid]) / denom
                temp2 <- sub.resp.mean[c(-mid,-mid - 1)]
                temp2 <- c(temp2, replace1)
            }
            
            ### collect the transformed variables into the vector ###
            
            temp <- c(temp, temp2)
        }
        
        ### calculate the absolute deviation from center ###
        
        resp.mean <- abs(temp)
        
        ### multiply the correction factor for the "zero.correction" option ###
        
        if (correction.method == "zero.correction")
        {
            correction <- sqrt((ngroup - 1) / ngroup)
            correction <- correction[-endpos]
            resp.mean <- resp.mean * correction
        }
    }
    
    orig.resp.mean <- resp.mean
    
    ### calculate p-values of the component test statistics ###
    
    for (i in 1:m)
    {
        resp.mean <- orig.resp.mean
        
        if (correction.method == "zero.correction" ||
            correction.method == "zero.removal")
        {
            n1 <- n[1:i] - 1
            j <- i + 1
            n2 <- n[j] - 1
        }
        
        else
        {
            n1 <- n[1:i]
            j <- i + 1
            n2 <- n[j]
        }
        sum1 <- sum(n1)
        sum2 <- sum(n2)
        sum3 <- sum1 + sum2
        subgroup <- c(rep(1, sum1), rep(2, sum2))
        sub.resp.mean <- resp.mean[1:sum3]
        mu <- mean(sub.resp.mean)
        z <- as.vector(sub.resp.mean - mu)
        d <- subgroup
        t.statistic <- summary(lm(z ~ d))$coefficients[2, 3]
        df <- summary(lm(z ~ d))$df[2]
        
        if (correlation.method == "pearson")
        {
            ### calculate the correlation between d and z ###
            
            correlation <- cor(d, z, method = "pearson")
            
            if (tail == "left")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(left-tailed with Pearson correlation coefficient)")
                p.value <-
                    pt(t.statistic, df, lower.tail = TRUE)
                log.p.value <- pt(t.statistic,
                                  df,
                                  lower.tail = TRUE,
                                  log.p = TRUE)
                log.q.value <- pt(t.statistic,
                                  df,
                                  lower.tail = FALSE,
                                  log.p = TRUE)
            }
            
            else if (tail == "right")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(right-tailed with Pearson correlation coefficient)")
                p.value <-
                    pt(t.statistic, df, lower.tail = FALSE)
                log.p.value <- pt(t.statistic,
                                  df,
                                  lower.tail = FALSE,
                                  log.p = TRUE)
                log.q.value <- pt(t.statistic,
                                  df,
                                  lower.tail = TRUE,
                                  log.p = TRUE)
            }
            
            else
            {
                tail = "both"
                
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(two-tailed with Pearson correlation coefficient)")
                p.value <-
                    pt(t.statistic, df, lower.tail = TRUE)
                log.p.value <- pt(t.statistic,
                                  df,
                                  lower.tail = TRUE,
                                  log.p = TRUE)
                log.q.value <- pt(t.statistic,
                                  df,
                                  lower.tail = FALSE,
                                  log.p = TRUE)
                
                if (p.value >= 0.5)
                {
                    p.value <- pt(t.statistic, df, lower.tail = FALSE)
                    log.p.value <- pt(t.statistic,
                                      df,
                                      lower.tail = FALSE,
                                      log.p = TRUE)
                    log.q.value <- pt(t.statistic,
                                      df,
                                      lower.tail = TRUE,
                                      log.p = TRUE)
                }
                
                p.value <- p.value * 2
                log.p.value <- log.p.value + log(2)
                log.q.value <- log(1 - p.value)
            }
        }
        
        else if (correlation.method == "kendall")
        {
            ### calculate the correlation between d and z ###
            
            correlation <- cor(d, z, method = "kendall")
            
            if (tail == "left")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(left-tailed with Kendall correlation coefficient)")
                p.value.temp <- Kendall(d, z)$sl
                
                if (correlation < 0)
                {
                    p.value <- p.value.temp / 2
                }
                
                else
                {
                    p.value <- 1 - p.value.temp / 2
                }
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
            
            if (tail == "right")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(right-tailed with Kendall correlation coefficient)")
                p.value.temp <- Kendall(d, z)$sl
                
                if (correlation > 0)
                {
                    p.value <- p.value.temp / 2
                }
                
                else
                {
                    p.value <- 1 - p.value.temp / 2
                }
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
            
            if (tail == "both")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(two-tailed with Kendall correlation coefficient)")
                p.value <- Kendall(d, z)$sl
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
        }
        
        else
        {
            ### calculate the correlation between d and z ###
            
            correlation <- cor(d, z, method = "spearman")
            
            if (tail == "left")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(left-tailed with Spearman correlation coefficient)")
                p.value.temp <-
                    cor.test(d, z, method = "spearman")$p.value
                if (correlation < 0)
                {
                    p.value <- p.value.temp / 2
                }
                
                else
                {
                    p.value <- 1 - p.value.temp / 2
                }
                
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
            
            if (tail == "right")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(right-tailed with Spearman correlation coefficient)")
                p.value.temp <-
                    cor.test(d, z, method = "spearman")$p.value
                if (correlation > 0)
                {
                    p.value <- p.value.temp / 2
                }
                
                else
                {
                    p.value <- 1 - p.value.temp / 2
                }
                
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
            
            if (tail == "both")
            {
                if (i == 1)
                    METHOD <-
                        paste(METHOD,
                              "(two-tailed with Spearman correlation coefficient)")
                p.value <-
                    cor.test(d, z, method = "spearman")$p.value
                q.value <- 1 - p.value
                log.p.value <- log(p.value)
                log.q.value <- log(q.value)
            }
        }
        
        ### store the log-transformed p-values (for computational purposes) ###
        
        p[i] <- log.p.value
        q[i] <- log.q.value
        non.bootstrap.p[i] <- log.p.value
        non.bootstrap.q[i] <- log.q.value
        
        ### perform bootstrapping ###
        
        if (bootstrap == TRUE)
        {
            if (i == 1)
                METHOD = paste("bootstrap", METHOD)
            
            ### step 2 of Lim and Loh (1996): initialize variables ###
            
            R <- 0
            n.trim <- n[1:j]
            y.trim <- y[1:sum(n.trim)]
            N.trim <- length(y.trim)
            group.trim <- group[1:sum(n.trim)]
            ngroup.trim <- ngroup[1:sum(n.trim)]
            orig.i <- i
            
            ### step 3 of Lim and Loh (1996): calculate the fractional trimmed mean ###
            
            frac.trim.alpha = 0.2
            b.trimmed.mean <- function(y)
            {
                nn <- length(y)
                wt <- rep(0, nn)
                y2 <- y[order(y)]
                lower <- ceiling(nn * frac.trim.alpha) + 1
                upper <- floor(nn * (1 - frac.trim.alpha))
                
                if (lower > upper)
                    stop("frac.trim.alpha value is too large")
                
                m <- upper - lower + 1
                frac <- (nn * (1 - 2 * frac.trim.alpha) - m) / 2
                wt[lower - 1] <- frac
                wt[upper + 1] <- frac
                wt[lower:upper] <- 1
                return(weighted.mean(y2, wt))
            }
            
            b.trim.means <-
                tapply(y.trim, group.trim, b.trimmed.mean)
            rm <- y.trim - b.trim.means[group.trim]
            
            ### step 7 of Lim and Loh (1996): enter a loop ###
            j <- 1
            
            for (j in 1:num.bootstrap)
            {
                orig.j <- j
                
                ### step 4 of Lim and Loh (1996): obtain a bootstrap sample ###
                
                sam <- sample(rm, replace = TRUE)
                boot.sample <- sam
                
                ### step 5 of Lim and Loh (1996): smooth the variables if n_i < 10 for at least one sample size ###
                
                if (min(n.trim) < 10)
                {
                    U <- runif(1) - 0.5
                    means <- tapply(y.trim, group.trim, mean)
                    v.trim <-
                        sqrt(sum((y.trim - means[group.trim]) ^ 2) / N.trim)
                    boot.sample <-
                        ((12 / 13) ^ (0.5)) * (sam + v.trim * U)
                }
                
                ### step 6 of Lim and Loh (1996): compute the bootstrap statistic, and increment R to R + 1 if necessary ###
                
                if (location == "mean") {
                    boot.means <- tapply(boot.sample, group.trim, mean)
                } else if (location == "median") {
                    boot.means <- tapply(boot.sample, group.trim, median)
                } else {
                    location <- "trim.mean"
                    boot.means <- tapply(boot.sample, group.trim, mean, trim = trim.alpha)
                }
                
                ### calculate bootstrap statistic ###
                z <- boot.sample - boot.means[group.trim]
                
                ### multiply the correction factor to each observation if "correction.factor" is chosen ###
                
                if (correction.method == "correction.factor")
                {
                    correction.trim <- sqrt(ngroup.trim / (ngroup.trim - 1))
                    z <- z * correction.trim
                }
                
                resp.boot.mean <- abs(z)
                
                ### perform correction techniques for "zero.removal" (Hines and Hines, 2000) ###
                ### or "zero.correction" (Noguchi and Gel, 2009).                            ###
                
                if (correction.method == "zero.correction" ||
                    correction.method == "zero.removal")
                {
                    ### set up variables for calculating the deviation from center ###
                    
                    resp.mean <- z
                    k <- length(n.trim)
                    temp <- double()
                    endpos <- double()
                    startpos <- double()
                    
                    ### calculate the absolute deviation from mean and remove zeros ###
                    
                    for (i in 1:k)
                    {
                        group.size <- n.trim[i]
                        j <- i - 1
                        
                        ### calculate the starting and ending index of each group ###
                        
                        if (i == 1)
                            start <- 1
                        else
                            start <- sum(n.trim[1:j]) + 1
                        startpos <- c(startpos, start)
                        end <- sum(n.trim[1:i])
                        endpos <- c(endpos, end)
                        
                        ### extract the deviation from center for the ith group ###
                        
                        sub.resp.mean <- resp.mean[start:end]
                        sub.resp.mean <-
                            sub.resp.mean[order(sub.resp.mean)]
                        
                        ### remove structural zero for the odd-sized group ###
                        
                        if (group.size %% 2 == 1)
                        {
                            mid <- (group.size + 1) / 2
                            temp2 <- sub.resp.mean[-mid]
                        }
                        ### remove structural zero for the even-sized group ###
                        
                        if (group.size %% 2 == 0)
                        {
                            mid <- group.size / 2
                            
                            ### set up the denominator value for the transformation ###
                            
                            ### set 1 for the "zero.correction" option ###
                            
                            denom <- 1
                            
                            ### set sqrt(2) for the "zero.removal" option ###
                            
                            if (correction.method == "zero.removal")
                                denom <- sqrt(2)
                            
                            ### perform the orthogonal transformation ###
                            
                            replace1 <-
                                (sub.resp.mean[mid + 1] - sub.resp.mean[mid]) / denom
                            temp2 <-
                                sub.resp.mean[c(-mid,-mid - 1)]
                            temp2 <- c(temp2, replace1)
                        }
                        
                        ### collect the transformed variables into the vector ###
                        
                        temp <- c(temp, temp2)
                    }
                    
                    ### calculate the absolute deviation from center ###
                    
                    resp.boot.mean <- abs(temp)
                    
                    ### multiply the correction factor for the "zero.correction" option ###
                    
                    if (correction.method == "zero.correction")
                    {
                        correction <- sqrt((ngroup.trim - 1) / ngroup.trim)
                        correction <- correction[-endpos]
                        resp.boot.mean <-
                            resp.boot.mean * correction
                    }
                }
                
                ### calculate the bootstrap statistic ###
                
                i <- orig.i
                
                if (correction.method == "zero.correction" ||
                    correction.method == "zero.removal")
                {
                    n1.trim <- n.trim[1:i] - 1
                    j <- i + 1
                    n2.trim <- n.trim[j] - 1
                } else {
                    n1.trim <- n.trim[1:i]
                    j <- i + 1
                    n2.trim <- n.trim[j]
                }
                
                sum1 <- sum(n1.trim)
                sum2 <- sum(n2.trim)
                sum3 <- sum1 + sum2
                subgroup <- c(rep(1, sum1), rep(2, sum2))
                sub.resp.mean <- resp.boot.mean[1:sum3]
                mu <- mean(sub.resp.mean)
                z <- as.vector(sub.resp.mean - mu)
                d <- subgroup
                correlation2 <-
                    cor(d, z, method = correlation.method)
                
                if (tail == "right")
                {
                    if (correlation2 > correlation)
                        R <- R + 1
                }
                
                else if (tail == "left")
                {
                    if (correlation2 < correlation)
                        R <- R + 1
                }
                
                else
                {
                    tail = "both"
                    
                    if (abs(correlation2) > abs(correlation))
                        R <- R + 1
                }
            }
            
            ### step 8 of Lim and Loh (1996): calculate the bootstrap p-value ###
            
            p.value <- R / num.bootstrap
            
            i <- orig.i
            
            p[i] <- log(p.value)
            q[i] <- log(1 - p.value)
        }
    }
    
    ### combine p-values using four p-value combining methods (T,F,N,L) ###
    
    psiT <- exp(min(p))
    psiF <- -2 * sum(p)
    psiN <- -sum(qnorm(p, lower.tail = TRUE, log.p = TRUE))
    A <- pi ^ 2 * m * (5 * m + 2) / (15 * m + 12)
    psiL <- -A ^ (-1 / 2) * sum(p - q)
    psiT.pvalue <- 1 - (1 - psiT) ^ m
    psiF.pvalue <- pchisq(psiF, df = 2 * m, lower.tail = FALSE)
    psiN.pvalue <- pnorm(psiN,
                         mean = 0,
                         sd = sqrt(m),
                         lower.tail = FALSE)
    psiL.pvalue <- pt(psiL, df = 5 * m + 4, lower.tail = FALSE)
    
    non.bootstrap.psiT <- exp(min(non.bootstrap.p))
    non.bootstrap.psiF <- -2 * sum(non.bootstrap.p)
    non.bootstrap.psiN <-
        -sum(qnorm(non.bootstrap.p, lower.tail = TRUE, log.p = TRUE))
    non.bootstrap.psiL <-
        -A ^ (-1 / 2) * sum(non.bootstrap.p - non.bootstrap.q)
    non.bootstrap.psiT.pvalue <-
        1 - (1 - non.bootstrap.psiT) ^ m
    non.bootstrap.psiF.pvalue <-
        pchisq(non.bootstrap.psiF,
               df = 2 * m,
               lower.tail = FALSE)
    non.bootstrap.psiN.pvalue <-
        pnorm(
            non.bootstrap.psiN,
            mean = 0,
            sd = sqrt(m),
            lower.tail = FALSE
        )
    non.bootstrap.psiL.pvalue <-
        pt(non.bootstrap.psiL,
           df = 5 * m + 4,
           lower.tail = FALSE)
    
    ### display output ###
    
    PSIT = psiT
    PSIF = psiF
    PSIN = psiN
    PSIL = psiL
    
    names(PSIT) = "Test Statistic (T)"
    names(PSIF) = "Test Statistic (F)"
    names(PSIN) = "Test Statistic (N)"
    names(PSIL) = "Test Statistic (L)"
    
    list(
        T = structure(
            list(
                statistic = PSIT,
                p.value = psiT.pvalue,
                data.name = DNAME,
                non.bootstrap.statistic = non.bootstrap.psiT,
                non.bootstrap.p.value = non.bootstrap.psiT.pvalue,
                method = METHOD
            ),
            class = "htest"
        ),
        F = structure(
            list(
                statistic = PSIF,
                p.value = psiF.pvalue,
                data.name = DNAME,
                non.bootstrap.statistic = non.bootstrap.psiF,
                non.bootstrap.p.value = non.bootstrap.psiF.pvalue,
                method = METHOD
            ),
            class = "htest"
        ),
        N = structure(
            list(
                statistic = PSIN,
                p.value = psiN.pvalue,
                data.name = DNAME,
                non.bootstrap.statistic = non.bootstrap.psiN,
                non.bootstrap.p.value = non.bootstrap.psiN.pvalue,
                method = METHOD
            ),
            class = "htest"
        ),
        L = structure(
            list(
                statistic = PSIL,
                p.value = psiL.pvalue,
                data.name = DNAME,
                non.bootstrap.statistic = non.bootstrap.psiL,
                non.bootstrap.p.value = non.bootstrap.psiL.pvalue,
                method = METHOD
            ),
            class = "htest"
        )
    )
}
