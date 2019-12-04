library(tidyverse)
library(here)

###########################################################################
### create Dataframe of BEDs ##############################################
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  Input - merged beds generated with mergeBEDs.sh 
###########################################################################

## create vector of file paths to be used for dataframes
all_38_beds <- list.files(here("data/stratifications/GRCh38/merged/"),
                          full.names = TRUE)
all_37_beds <- list.files(here("data/stratifications/GRCh37/merged/"),
                          full.names = TRUE)

## Create dataframe of file data using list files
all_38_beds_set <- set_names(all_38_beds)
all_37_beds_set <- set_names(all_37_beds)

all_38_beds_df <- map_dfr(.x= all_38_beds_set, read_tsv,
                          col_types = "cdd",
                          col_names= c("chrom","start","end"),
                          .id = "stratification_filename") %>% 
  #change stratification column to just the base name of file, removing path
  mutate(stratification_filename = basename(stratification_filename)) 

all_37_beds_df <- map_dfr(.x= all_37_beds_set, read_tsv,
                          col_names= c("chrom","start","end"),
                          col_types = "cdd",
                          .id = "stratification_filename") %>% 
  #change stratification column to just the base name of file, removing path
  mutate(stratification_filename = basename(stratification_filename))

## Jenny note-to-self ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## `col_types` need to be specified because `readr` tries to use the first 1000
## rows to infer the datatype. GRCh37 chromsome go from 22 -> X therefore
## datatype changes.  
##
## `col_types = "cdd"` set the data types for the columns to character, 
## double, double.
##
## Issue with refseq_union_cds.sort.bed.gz, this file turns out was empty
## and was throwing an error, 
##     """"
##     Error in read_tokens_(data, tokenizer, col_specs, col_names,
##     locale_:Cannot read file
##     /Users/jmcdani/Documents/GiaB/Benchmarking/GRCh38_stratification_validation/data/stratifications/GRCh37/merged/foo:
##     Invalid argument
##    """"
## Nate helped identify the problem file using safe_tsv <- safely(read_tsv) and
## transpose functions. 
##This was only an issue with the GRCh37 refseq stratification.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Joining Dataframes using stratfile_summary

## join dataframes of beds and stratification levels
stratfile_summary <- here("data/stratification_summary_spreadsheets",
                          "110619_37and38_stratification_types_levels_md5_longFormat.csv") %>% 
  read.csv(header = TRUE, na.strings = TRUE)

all_38_beds_df_with_levels <- left_join(all_38_beds_df, stratfile_summary)
all_37_beds_df_with_levels <- left_join(all_37_beds_df, stratfile_summary)

###########################################################################
## Chromosome Presence/ Absence Plots #####################################
###########################################################################
## vector of chromosomes we expect to have
chroms <- c(1:22, "X", "Y")
expected_chroms <- factor(c(paste0("chr", chroms), chroms), 
                          levels = c(paste0("chr", chroms), chroms))
                          


## dataframes with Expected and Unexpected chroms for plotting
expected_chroms_38 <- all_38_beds_df_with_levels %>%
  filter(chrom %in% expected_chroms)
unexpected_chroms_38 <- all_38_beds_df_with_levels %>%
  filter(!chrom %in% expected_chroms)

expected_chroms_37 <- all_37_beds_df_with_levels %>%
  filter(chrom %in% expected_chroms)
unexpected_chroms_37 <- all_37_beds_df_with_levels %>%
  filter(!chrom %in% expected_chroms)

#plot expected chroms
expected_chroms_38 %>%
  select(chrom, stratification_type, stratification_level) %>%
  distinct() %>%
  ggplot(aes(x = stratification_level, y = chrom)) +
  geom_point(aes(col = stratification_type), show.legend = FALSE) +
  labs(x = "Stratification", y = "Chrom",
       title = "Expected Chromosome Presence 38 Stratifications") +
  theme(axis.text.x = element_text(angle = 90))

expected_chroms_37 %>%
  select(chrom, stratification_type, stratification_level) %>%
  distinct() %>%
  ggplot(aes(x = stratification_level, y = chrom)) +
  geom_point(aes(col = stratification_type), show.legend = FALSE) +
  labs(x = "stratification", y = "chrom",
       title = "Chromosome Presence 37 Stratifications") +
  theme(axis.text.x = element_text(angle = 90))

# plot UNexpected chroms
unexpected_chroms_38 %>%
  select(chrom, stratification_type, stratification_level) %>%
  distinct() %>% 
  ggplot(aes(x = stratification_level, y = chrom)) +
  geom_point(aes(col = stratification_type), show.legend = FALSE) +
  labs(x = "stratification", y = "chrom",
       title = "UNexpected Chromosome Prescence 38 stratifications") +
  theme(axis.text.x = element_text(angle = 90))

unexpected_chroms_37 %>%
           select(chrom, stratification_type, stratification_level) %>%
  distinct() %>% 
  ggplot(aes(x= stratification_level,y = chrom)) +
  geom_point(aes(col=stratification_type), show.legend = FALSE) +
  labs(x="stratification", y="chrom",
       title="UNexpected Chromosome Prescence 37 stratifications") +
  theme(axis.text.x = element_text(angle = 90))

###########################################################################
## Get Chromosome Sizes ###################################################
###########################################################################

## Jenny note-to-self ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Nate's Function from v3.3.2 paper github --> R --> functions
## https://github.com/nate-d-olson/giab-integrationV3.3.2-paper/blob/master/R/functions.R
## Get non-N base genome and chromosome sizes,
## requires BSgenome formatted object
## to inspect genome object create a genome obj using genome <- get_genome("genome name")
## you can see obj by using things like seqnames(genome) to look at chroms, to look
## at particular elements use genome[[element #]], to see what is in a particular
## element seqnames(genome)[["element #]]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_genome <- function(genome){
  if (genome == "hs37d5"){
    require(BSgenome.Hsapiens.1000genomes.hs37d5)
    return(BSgenome.Hsapiens.1000genomes.hs37d5)
  } else if (genome == "GRCh38") {
    require(BSgenome.Hsapiens.NCBI.GRCh38)
    return(BSgenome.Hsapiens.NCBI.GRCh38)
  } else {
    stop("Genome not `hs37d5` or `GRCh38`")
  }
}

get_chrom_sizes <- function(genome){

  genome_obj <- get_genome(genome)

  get_alpha_freq <- function(i){
    genome_obj[[i]] %>%
      alphabetFrequency() %>%
      data.frame()
  }

  ## Adjusted list to accomodate 24 elements to include X and Y.  
  ## Original code only used chroms 1:22, `alpha_freq_df <- as.list(1:22)`
  c_list <- as.list(as.character(1:24))
  c_list[23:24] <- c("X","Y")
  alpha_freq_df <- c_list %>%
    map_dfc(get_alpha_freq)

  colnames(alpha_freq_df) <- paste0("chr",1:24)

  alpha_freq_df <- alpha_freq_df %>%
    ## Removing bases not included in counts and non-standard bases
    filter(chr1 > 100) %>% ## Not sure what this is doing here....
    add_column(base = c("A","C","G","T","N"))

  chromosome_lengths <- alpha_freq_df %>%
    tidyr::gather(key = "chrom", value = "nbases", -base) %>%
    group_by(chrom) %>%
    mutate(base_type = if_else(base == "N", "N", "non_N")) %>%
    group_by(chrom, base_type) %>%
    summarise(n_bases = sum(nbases)) %>%
    tidyr::spread(base_type, n_bases) %>%
    mutate(len = N + non_N) %>%
    dplyr::select(-N)

  ## data frame with total length and number of non-N bases
  data_frame(chrom = "genome",
             non_N = sum(chromosome_lengths$non_N),
             len = sum(chromosome_lengths$len)) %>%
    bind_rows(chromosome_lengths)
}

get_genome("hs37d5")
chromsizes_37 <- get_chrom_sizes("hs37d5")
## change chr23 and chr24 to chrX and chrY
chromsizes_37[17,1] <- c("chrX")
chromsizes_37[18,1] <- c("chrY")

## BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38") - only do first time
get_genome("GRCh38")
chromsizes_38 <- get_chrom_sizes("GRCh38")
## change chr23 and chr24 to chrX and chrY
chromsizes_38[17,1] <- c("chrX")
chromsizes_38[18,1] <- c("chrY")


###########################################################################
## Create Dataframes for plotting chromosome sizes ########################
###########################################################################

## Add column to table of 37 and 38  to calculate region sizes
expected_chroms_37 <- mutate(expected_chroms_37, region_size = end - start)
expected_chroms_38 <- mutate(expected_chroms_38, region_size = end - start)

## Jenny note-to-self ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Needed to find way to group and sum by filename and group for each filename and
## chrom combination sum region_size column which was helpful to find which row
## chrom1 ended and chrom2 began. Used to look at summing and how it was working.
## which(expected_chroms_37$stratification_filename ==
## "AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz" & expected_chroms_37$chrom == "2")
## Nate helped with below using group_by and summarise functions
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Sum of all regions by chromosome
chrom_region_sizes_37stratifications <- expected_chroms_37 %>%
  group_by(stratification_filename, chrom) %>%
  summarise(total_region = sum(region_size))

chrom_region_sizes_38stratifications <- expected_chroms_38 %>%
  group_by(stratification_filename, chrom) %>%
  summarise(total_region = sum(region_size))


## append "chr" prefix to chroms in GRCh37 df to reduce number of fields when
## joined with 38. Found there were some with "chr" already appended therefore
## needed to remove all chr first then paste chr to ensure single chr prefix

chrom_region_sizes_37stratifications$chrom <-
  str_remove_all(chrom_region_sizes_37stratifications$chrom, "chr")
chrom_region_sizes_37stratifications$chrom <-
  paste0("chr", chrom_region_sizes_37stratifications$chrom)

## combine ref and strat region size dataframes
## total_region=sum of stratification chrom regions size (use for plots),
## Non_N=ref chrom size w/o Ns (use for plots), len=total length of chrom regions

stratANDref_37chromsizes <- left_join(chrom_region_sizes_37stratifications,
                                      chromsizes_37)
stratANDref_37chromsizes <- add_column(stratANDref_37chromsizes,
                                       "reference" = "GRCh37")

stratANDref_38chromsizes <- left_join(chrom_region_sizes_38stratifications,
                                      chromsizes_38)
stratANDref_38chromsizes <- add_column(stratANDref_38chromsizes,
                                       "reference" = "GRCh38")


## Create table with all stratification and reference chrom sizes and
## stratification types
full <- full_join(stratANDref_37chromsizes,stratANDref_38chromsizes)

summary_for_plots <- full %>% left_join(stratfile_summary) %>%
  rename(non_N = "non_N_ref", len = "len_ref",
         total_region = "total_region_stratification") %>%
  mutate(coverage = total_region_stratification/non_N_ref)

###########################################################################
## Coverage Plots #########################################################
###########################################################################
coverage_plot <- function(strat_type_df, strat_type){
  ggplot(strat_type_df, 
         aes(x=chrom, y=coverage, group = reference)) +
    geom_col(aes(fill = reference), position = "dodge") +
    facet_grid(scales = "free_y", stratification_level ~ .) +
    theme(axis.text.x = element_text(angle = 90), 
          legend.position = "top",
          strip.text.y = element_text(angle = 0)) +
    labs(title = strat_type)
}

summary_for_plots %>% 
  filter(stratification_type == "FunctionalRegions") %>% 
  coverage_plot(strat_type = "FunctionalRegions")

summary_for_plots %>% 
  filter(stratification_type == "GCcontent") %>% 
  coverage_plot(strat_type = "GCcontent")

summary_for_plots %>% 
  filter(stratification_type == "mappability") %>% 
  coverage_plot(strat_type = "mappability")

summary_for_plots %>% 
  filter(stratification_type == "SegmentalDuplications") %>% 
  coverage_plot(strat_type = "SegmentalDuplications")

summary_for_plots %>% 
  filter(stratification_type == "FunctionalTechnicallyDifficultRegions") %>% 
  coverage_plot(strat_type = "FunctionalTechnicallyDifficultRegions")

summary_for_plots %>% 
  filter(stratification_type == "LowComplexity") %>% 
  coverage_plot(strat_type = "LowComplexity")


GSpecific <-filter(summary_for_plots, stratification_type == "GenomeSpecific")
ggplot(GSpecific, aes(x=chrom, y=coverage, group = reference)) +
  geom_col(aes(fill = reference ), position = "dodge") +
  facet_grid(stratification_level ~ genome, scales = "free") +
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = "top",
        strip.text.y = element_text(angle = 0)) +
  labs(title = "GenomeSpecific")

