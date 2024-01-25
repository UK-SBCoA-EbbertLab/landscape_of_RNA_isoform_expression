sink('../../tables/GTEx_expression/proportion_test_statistic.txt')
prop.test(x=c(22522-14649, 12989-10660), n=c(22522, 12989))
sink()

