library(tidyverse)

p_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_preferential_filtered.tsv') %>%
	group_by(TXNAME) %>%
	group_split()

print(p_iso)

p_isoforms <- tibble(
		     TXNAME = character(),
		     type = character())

for (tab in p_iso) {
	tmp <- c()
	tmp <- append(tmp, tab$tis1)
	tmp <- append(tmp, tab$tis2)
	summarized <- as.data.frame(table(tmp))
	n_sig <- unique(tab$n_sig)
	tmp8 <- sum(summarized$Freq == 8)
	tmp7 <- sum(summarized$Freq == 7)
	if (any(summarized$Freq == 8)) {
		p_isoforms <- p_isoforms %>%
			add_row(TXNAME = unique(tab$TXNAME),
		#		type = ifelse(n_sig == 8, 'one_tis', ifelse(n_sig == 15, 'two_tis', 'three_tis')))
				type = ifelse(tmp8 == 3 | (tmp8 == 1 & tmp7 == 2), 'three_tis', 
					       ifelse(tmp8 == 2, 'two_tis', 'one_tis')))
	} else if (sum(summarized$Freq == 7) == 2 & n_sig == 14) {
		p_isoforms <- p_isoforms %>%
			add_row(TXNAME = unique(tab$TXNAME),
				type = 'two_tis')
	} else if (sum(summarized$Freq == 6) == 3 & n_sig == 18) {
		p_isoforms <- p_isoforms %>%
			add_row(TXNAME = unique(tab$TXNAME),
				type = 'three_tis')
	}
}

p_iso <- read_tsv('../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_preferential_filtered.tsv') 

write_tsv(p_isoforms %>% left_join(p_iso) %>% distinct(), '../../tables/GTEx_expression_our_new_isoforms/deseq_output_all_new_isoforms_preferential_filtered_twice.tsv')
			
