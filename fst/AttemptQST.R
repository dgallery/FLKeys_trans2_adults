devtools::install_github('famuvie/breedR')

library(breedR)
?remlf90

#rows are genes, columns are individuals
#metadata that populations assigned to them 

#model gets fit to each gene... need a vector(?) of genes, how do we get this from a matrix?
#summed? Expression = overall expression 

#y is a gene expression vector across individuals... so is it a matrix?
#b and w are respectively random effects of populations and individuals within populations

#can calculate var in GE for a particular gene within a pop
#can calulate var in GE across population for a gene 

meta_gene = coldata
vsd.t = t(vsd)
vsd.df = as.data.frame(vsd)
vsd.df$geneID = rownames(vsd.df)
meta_gene$ID = rownames(meta_gene)

by_pop = pivot_longer(vsd.df, cols = !geneID, names_to = "ID", values_to = "expression") %>% left_join(meta_gene) %>%
  group_by(habitat, geneID) %>% summarise(exp_var = var(expression))

overall_near_off = pivot_longer(vsd.df, cols = !geneID, names_to = "ID", values_to = "expression") %>% left_join(meta_gene) %>%
  filter(habitat != "deep") %>%
  group_by(geneID) %>% summarise(exp_var = var(expression))

overall_near_deep = pivot_longer(vsd.df, cols = !geneID, names_to = "ID", values_to = "expression") %>% left_join(meta_gene) %>%
  filter(habitat != "offshore") %>%
  group_by(geneID) %>% summarise(exp_var = var(expression))

overall_off_deep = pivot_longer(vsd.df, cols = !geneID, names_to = "ID", values_to = "expression") %>% left_join(meta_gene) %>%
  filter(habitat != "nearshore") %>%
  group_by(geneID) %>% summarise(exp_var = var(expression))

####### ATTEMPT TO CALCULATE BY GENE Q_ST,
#assuming that between pop variance is just overall variance in expression
qst = function(b_var, w_var){
  qst = b_var^2/(b_var^2 + 2*(w_var^2))
  return(qst)
}

#big question is how to calculate sigma_b and sigma_w. It feels like you need to account for
#the genetic variation 

#are we assuming that within population variances are the same for each gene????
#is that why it's times 2 

