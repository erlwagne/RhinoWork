# Function to automate pairwise comparison of a list of loo objects
# Value is a list whose i-th element is the pairwise comparison returned by 
#   loo::compare(loos[i], loos[i+1])

loo_compair <- function(loos) {
  compr <- loo::loo_compare(loos)
  compair <- vector("list", length(loos) - 1) 
  for(i in 1:length(compair))
  {
    names(compair)[i] <- paste(row.names(compr)[c(i+1,i)], collapse = " vs. ")
    compair[[i]] <- 2*loo::compare(x = loos[row.names(compr)[c(i+1,i)]])
    names(compair[[i]])[1] <- "looic_diff"
  }
  return(compair)
}