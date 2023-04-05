detect_family <- function(Z, family = NULL){
  if( inherits(family, "family") ||  (!is.null(family) && family != "auto") ) return(family)
  if(all(Z %in% c(0,1))) return("binomial")
  if(all(Z >=0)) return("poisson")
  return("gaussian")
}
