sink('../../tables/GTEx_expression/proportion_test_statistic.txt')
results = prop.test(x=c(22522-14649, 12989-10660), n=c(22522, 12989))
print(results)
results$p.value

# chi-square test of independence 
table <- matrix(c(14649, 22522-14649, 
		  10660, 12989-10660),
		  nrow=2,
		  dimnames=list(c('Cerebellar Hemisphere', 'Left Ventricle'), c('Protein-coding', 'Non protein-coding')))
chi<-chisq.test(table)
print(chi)
chi$p.value
fish <-fisher.test(table)
print(fish)
fish$p.value
sink()

