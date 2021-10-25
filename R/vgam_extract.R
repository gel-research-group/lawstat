#These function were taken from VGAM (v. 1.1-2) because the package was going to be archived.
# Thomas W. Yee (2019). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-2. URL
# https://CRAN.R-project.org/package=VGAM


VGAM_plaplace = function (q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) 
{
    zedd <- (q - location)/scale
    if (!is.logical(lower.tail) || length(lower.tail) != 1) 
        stop("bad input for argument 'lower.tail'")
    if (!is.logical(log.p) || length(log.p) != 1) 
        stop("bad input for argument 'log.p'")
    L <- max(length(q), length(location), length(scale))
    if (length(q) != L) 
        q <- rep_len(q, L)
    if (length(location) != L) 
        location <- rep_len(location, L)
    if (length(scale) != L) 
        scale <- rep_len(scale, L)
    if (lower.tail) {
        if (log.p) {
            ans <- ifelse(q < location, log(0.5) + zedd, log1p(-0.5 * exp(-zedd)))
        }
        else {
            ans <- ifelse(q < location, 0.5 * exp(zedd), 1 - 0.5 * exp(-zedd))
        }
    }
    else {
        if (log.p) {
            ans <- ifelse(q < location, log1p(-0.5 * exp(zedd)), 
                          log(0.5) - zedd)
        }
        else {
            ans <- ifelse(q < location, 1 - 0.5 * exp(zedd), 
                          0.5 * exp(-zedd))
        }
    }
    ans[scale <= 0] <- NaN
    ans
}