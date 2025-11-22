bern_prob <- function(rho, p, q) {
p00 = rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
p10 = 1 - q - p00
p01 = 1 - p - p00
p11 = p + q + p00 -  1
p00 = max(0, p00); p10 = max(0, p10); p01 = max(0, p01); p11 = max(0, p11)
return(c(p00, p01, p10, p11))
}

