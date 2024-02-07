library(tidyverse)

cpm <- read_tsv('../../tables/GTEx_expression/GTEx_isoforms_in_tissues_passing_med_CPM_gt_1.tsv')

median(cpm$Liver, na.rm=TRUE)
mean(cpm$Liver, na.rm=TRUE)

median(cpm$Lung, na.rm=TRUE)
mean(cpm$Lung, na.rm=TRUE)

median(cpm %>% pull(`Brain - Cerebellar Hemisphere`), na.rm=TRUE)
mean(cpm %>% pull(`Brain - Cerebellar Hemisphere`), na.rm=TRUE)

median(cpm %>% pull(`Heart - Left Ventricle`), na.rm=TRUE)
mean(cpm %>% pull(`Heart - Left Ventricle`), na.rm=TRUE)

median(cpm %>% pull(`Muscle - Skeletal`), na.rm=TRUE)
mean(cpm %>% pull(`Muscle - Skeletal`), na.rm=TRUE)
