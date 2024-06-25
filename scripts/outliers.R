# From Brenna R. Forester at: https://popgen.nescent.org/2018-03-27_RDA_GEA.html
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}