library(tidyverse)

#' Just a quick try to find perfect co-occurrence and co-exclusion. 
#' Easy to understand, implement and calculate n.


#' Find perfect co-occurrence
d <- tibble(a = c(1, 0, 1, 0),
            b = c(1, 0, 1, 0)) %>%
            filter(a > 0 | b > 0) 

if( nrow(d) == 0 ) {
  co <- FALSE
} else {
  co <- all( ( (d$a > 0) + (d$b > 0) ) == 2 )
}
co


#' Find perfect co-exclusion
d <- tibble(a = c(1, 0, 0, 0),
            b = c(0, 0, 0, 1)) %>%
            filter(a > 0 | b > 0) 

if( nrow(d) == 0 ) {
  ex <- FALSE
} else if(all(d$a == 0) || all(d$b == 0) ) {
  ex <- FALSE
} else {
  ex <- all(xor(d$a, d$b))
}
ex
