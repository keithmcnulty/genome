library(readr)
library(dplyr)
library(gwascat)
library(stringr)
library(ggplot2)

updated_gwas_data <- as.data.frame(makeCurrentGwascat())

my_genome_file <- "genome_Keith_McNulty_v3_Full_20191031074721.txt" ## <-- enter path to your genome file (23andme or ancestry should work)

my_genome_data <- readr::read_tsv(my_genome_file, skip = 19, col_types = readr::cols(
  `# rsid` = col_character(),
  chromosome = col_character(),
  position = col_integer(),
  genotype = col_character())
)

my_genome_data <- my_genome_data %>% 
  dplyr::rename(rsid = `# rsid`)

output_data <- inner_join(my_genome_data, updated_gwas_data, by = c("rsid" = "SNPS"))

output_data$risk_allele_clean <- str_sub(output_data$STRONGEST.SNP.RISK.ALLELE, -1)
output_data$my_allele_1 <- str_sub(output_data$genotype, 1, 1)
output_data$my_allele_2 <- str_sub(output_data$genotype, 2, 2)
output_data$have_risk_allele_count <- if_else(output_data$my_allele_1 == output_data$risk_allele_clean, 1, 0) + if_else(output_data$my_allele_2 == output_data$risk_allele_clean, 1, 0)

output_data %>% 
  dplyr::select(rsid, DISEASE.TRAIT, risk_allele = risk_allele_clean, your_geneotype = genotype) %>% 
  dplyr::filter(str_detect(tolower(DISEASE.TRAIT), "omega") | str_detect(tolower(DISEASE.TRAIT), "fatty acid")) %>%
  dplyr::distinct()

trait_entry_count <- output_data %>% 
  dplyr::group_by(DISEASE.TRAIT) %>%
  dplyr::filter(have_risk_allele_count >= 1) %>%
  dplyr::summarise(count_of_entries = n())

ggplot(filter(trait_entry_count, grepl("MTAG", DISEASE.TRAIT)), aes(x = reorder(DISEASE.TRAIT, count_of_entries, sum), y = count_of_entries)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(title = "Risk allele entry count for MTAG studies", y = "Count of MTAG entries", x = "Trait")
