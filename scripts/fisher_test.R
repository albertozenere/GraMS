# Auxiliary function used by multiple function to execute Fisher's exact test


fisher_test <- function(list_1, background ,list_2){
  
  stopifnot(class(list_1)=="character")
  stopifnot(class(background)=="character")
  stopifnot(class(list_2)=="character")
  
  list_1 <- unique(list_1)
  background <- unique(background)
  list_2 <- unique(list_2)
  
  stopifnot(all(list_1 %in% background))
  stopifnot(all(list_2 %in% background))
  
  
  #contingency table
  a <- list_1 %>% intersect(list_2) %>% length()
  b <- list_1 %>% setdiff(list_2) %>% length()
  c <- list_2 %>% setdiff(list_1) %>% length()
  d <- length(background) - a - b- c
  
  cont_mat <- matrix(c(a,b,c,d), nrow=2, ncol=2, byrow = TRUE)
  
  #test
  test <- fisher.test(cont_mat, alternative="greater")
  
  #output result
  out <- data.frame(p_value = test$p.value, odds_ratio = test$estimate["odds ratio"], n_genes=a,
                    row.names = NULL)
  
  
  return(out)
  
  
}
