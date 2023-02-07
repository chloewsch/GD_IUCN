## Code accompanying Schmidt et al. (2023) Genetic diversity and IUCN Red List status. Conservation Biology
## 1. Data prep

# Libraries--------
library(tidyr)
library(dplyr)
library(rredlist)

# load redlist api key

# microsatellite data ------------
## Data from: Lawrence et al. 2019 Scientific Data, Geo-referenced population-specific microsatellite data across American continents, the MacroPopGen Database
## https://www.nature.com/articles/s41597-019-0024-7
mpg <- read.csv('Macropopgen.csv', h=T,
                na.strings = "NA")
#load('IUCN_results_08_Mar_22b.RData') # IUCN data as of March 8

## with rredlist ####
mpg$row_id <- c(1:nrow(mpg))

# Create column for IUCN search
mpg$to_IUCN <- gsub('_', ' ', mpg$G_s)

# Remove subspecies
mpg$to_IUCN <- word(mpg$to_IUCN, 1, 2, sep = " ") # keeps only first 2 words

# Edit species names
mpg$to_IUCN <- gsub('Vireo atricapillus', 'Vireo atricapilla', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Pantherophis obsoleta', 'Pantherophis obsoletus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Cychlura cychlura', 'Cyclura cyclura', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Xerospermosphilus polionotus', 'Peromyscus polionotus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Xerospermosphilus perotensis', 'Xerospermophilus perotensis', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Chlorospingus ophthalmicus', 'Chlorospingus flavopectus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Myrmeciza exsul', 'Poliocrania exsul', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Pipra filicaud', 'Pipra filicauda', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Volatinia jacarin', 'Volatinia jacarina', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Somateria fisher', 'Somateria fischeri', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Coturnicops noveboracencis', 'Coturnicops noveboracensis', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Vermivora pinus', 'Vermivora cyanoptera', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Hemignathus virens', 'Chlorodrepanis virens', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Serpentina serpentina', 'Chelydra serpentina', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Craugastor bransfordi', 'Craugastor bransfordii', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Oophaa pumilio', 'Oophaga pumilio', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Empidonax trailii', 'Empidonax traillii', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Empidonax trailliis', 'Empidonax traillii', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Charadriu melodus', 'Charadrius melodus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Procnias tricarunculata', 'Procnias tricarunculatus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Buarremon torquatus', 'Arremon torquatus', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Pipilo crissalis', 'Melozone crissalis', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Ammodramus nelsoni', 'Ammospiza nelsoni', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Hemignathus virens', 'Chlorodrepanis virens', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Perca flavesecens', 'Perca flavescens', mpg$to_IUCN)
mpg$to_IUCN <- gsub('Sander S.vitreuss', 'Sander vitreus', mpg$to_IUCN)
mpg$to_IUCN <- gsub(',', '', mpg$to_IUCN)

# Some species are only defined at Genus level; make NA:
no_spp <- mpg[grep('spp', mpg$G_s),]$row_id
no_sp <- mpg[grep(' sp$', mpg$to_IUCN),]$row_id
no_sps <- mpg_fish[grep(' sp\\.', mpg_fish$to_IUCN),]$row_id
no_S <- which(nchar(word(mpg_fish$to_IUCN, 2))<=2) # 2 rows where species == "S."

mpg[c(no_spp, no_sp, no_sps, no_S), 'to_IUCN'] <- NA

## Search ####
# Create list of species names to search
mpg_sp <- unique(mpg$to_IUCN)
mpg_sp <- na.omit(mpg_sp)
species_list <- as.list(mpg_sp)

# Search over species list
sp_iucn_list <- lapply(species_list, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Create dataframe
sp_iucn <- do.call(rbind, sp_iucn_list)

## Find species with no search results ####
no_res <- which(unlist(lapply(sp_iucn_list, length))==0)
sp_not_found <- species_list[no_res]

## Find synonyms and accepted names if possible ####
syn_list <- lapply(sp_not_found, function(x) rl_synonyms(x, key = IUCN_REDLIST_KEY)$result)
syn_df <- do.call(rbind, syn_list)

# See if any had no synonyms:
no_res_syn <- which(unlist(lapply(syn_list, length))==0)
syn_not_found <- unlist(sp_not_found[no_res_syn]) 

# Create new list of species to search without duplicates
syn_search <- unique(syn_df$accepted_name)

# Search:
syn_search_res <- lapply(syn_search, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Dataframe with synonym results:
syn_res_df <- do.call(rbind, syn_search_res)

## unified taxonomy ####
synonyms <- syn_df %>% ## (!) note syn_df has duplicate accepted names with multiple synonyms
  select(accepted_name, synonym) # synonym is the column that matches with macropopgen

# create a dataframe with accepted & synonym to merge all results:
syn_merge <- data.frame(accepted_name = sp_iucn$scientific_name,
                        synonym = sp_iucn$scientific_name)

synonyms <- rbind(synonyms, syn_merge)
names(synonyms)[1] <- "scientific_name"

## Clean specific ones:
# Scinax perpusillus matched to Ololygon perpusilla and Ololygon peixotoi; remove 2nd
grep('Scinax perpusillus', synonyms$synonym) # remove row 11
View(synonyms[c(grep('Scinax perpusillus', synonyms$synonym)),])
synonyms[synonyms$synonym == 'Scinax perpusillus', "scientific_name"] <- "Ololygon perpusilla"

# Cebus apella (reference R806) is the subspecies Cebus apella nigritus, AKA  Sapajus nigritus
grep('Cebus apella', synonyms$synonym)
View(synonyms[c(grep('Cebus apella', synonyms$synonym)),])
synonyms[synonyms$synonym == 'Cebus apella', "scientific_name"] <- "Sapajus nigritus"

## Remove duplicate rows:
synonyms <- distinct(synonyms)

## All results ####
all.res <- rbind(sp_iucn, syn_res_df)
all.res <- all.res %>% 
  select(taxonid, scientific_name, kingdom, phylum, class, order, family,
         genus, category, population_trend, marine_system, freshwater_system,
         terrestrial_system) %>% # keep only certain columns
  full_join(., synonyms, by = 'scientific_name') # add synonym column to match with macropopgen

## Merge with macropopgen ####
mpg_iucn <- merge(mpg, all.res, by.x = "to_IUCN", by.y = "synonym", all = TRUE, 
                  incomparables = NA)

mpg_iucn2 <- mpg_iucn %>% 
  drop_na(category) # drop species with no category

## Check instances where synonym mapped to multiple accepted names
check <- distinct(mpg_iucn2, to_IUCN, category)
check[c(which(duplicated(check$to_IUCN))),]

View(mpg_iucn2[c(which(is.na(mpg_iucn2$to_IUCN))),])

## Remove entries where to_IUCN is NA
mpg_iucn2 <- mpg_iucn %>% 
  drop_na(to_IUCN) %>% 
  filter(to_IUCN != 'Chelonoidis nigra')

## Genetic diversity summary ####
mpg_iucn3 <- mpg_iucn2

# Clean up categories
mpg_iucn3$rlcat <- gsub("LR/lc|LR/cd|LR/nt", "LC", mpg_iucn3$category)

mpg_iucn3 <- mpg_iucn3 %>%
  distinct(PopId, .keep_all = TRUE) %>%  # remove duplicate populations
  drop_na(category) %>% 
  group_by(scientific_name) %>% 
  mutate(meanGD = mean(He, na.rm = TRUE), # mean per species
         num_indiv = sum(n)) %>%
  ungroup() %>% 
  distinct(TaxaClass, scientific_name, meanGD, rlcat, population_trend, .keep_all = T) %>% 
  mutate(rlcat = factor(rlcat, levels = c('CR',
                                          'EN',
                                          'VU',
                                          'NT',
                                          'LC',
                                          'DD')))

mpg_iucn4 <- mpg_iucn3 %>%
  distinct(scientific_name, .keep_all = TRUE) %>% 
  drop_na(meanGD, rlcat) %>% # remove rows with missing genetic diversity or RL category
  filter(rlcat != 'DD') %>% # remove data deficient species
  mutate(IUCN_fact = as.factor(case_when(rlcat == 'CR' ~ '5',
                                         rlcat == 'EN' ~ '4',
                                         rlcat == 'VU' ~ '3',
                                         rlcat == 'NT' ~ '2',
                                         rlcat == 'LC' ~ '1')),
         IUCN_binary = ifelse(rlcat %in% c('CR', 'EN', 'VU'), 
                              1, 0))

## file -----
## Remove '\n' in cols: Common, G_s, and Location to write file
mpg_iucn4$Location <- gsub('\n', '', mpg_iucn4$Location, fixed = TRUE)
mpg_iucn4$G_s <- gsub('\n', '', mpg_iucn4$G_s, fixed = TRUE)
mpg_iucn4$Common <- gsub('\n', '', mpg_iucn4$Common, fixed = TRUE)

#write.table(mpg_iucn4, "microsat_data.txt", row.names = F, quote = F, sep = '\t')

# mitochondrial DNA data ------------
## Data from: Canteri et al. 2021 Ecography, IUCN Red List protects avian genetic diversity
## https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05895
mtdna <- read.csv("Canteri_mtdna.csv", header = T)

## Search ####
mtdna$SP <- gsub('_', ' ', mtdna$SP)
mtdna_sp <- unique(mtdna$SP)
species_listmt <- as.list(mtdna_sp)

# Search over species list
sp_iucn_list_mt <- lapply(species_listmt, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Create dataframe
sp_iucnmt <- do.call(rbind, sp_iucn_list_mt)

## Find species with no search results ####
no_resmt <- which(unlist(lapply(sp_iucn_list_mt, length))==0)
sp_not_foundmt <- species_listmt[no_resmt]

## Find synonyms and accepted names if possible ####
syn_listmt <- lapply(sp_not_foundmt, function(x) rl_synonyms(x, key = IUCN_REDLIST_KEY)$result)
syn_dfmt <- do.call(rbind, syn_listmt)

# Create new list of species to search without duplicates
syn_searchmt <- unique(syn_dfmt$accepted_name)

# Search:
syn_search_resmt <- lapply(syn_searchmt, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Dataframe with synonym results:
syn_res_dfmt <- do.call(rbind, syn_search_resmt)

## unified taxonomy ####
synonyms_mt <- syn_dfmt %>% ## (!) note syn_df may have duplicate accepted names with multiple synonyms
  dplyr::select(accepted_name, synonym) # synonym is the column that matches with original data

# create a dataframe with accepted & synonym to merge all results:
syn_mergemt <- data.frame(accepted_name = sp_iucnmt$scientific_name,
                          synonym = sp_iucnmt$scientific_name)

synonyms_mt <- rbind(synonyms_mt, syn_mergemt)
names(synonyms)[1] <- "scientific_name"

nrow(synonyms_mt)
length(unique(synonyms_mt$accepted_name))
length(unique(synonyms_mt$synonym))

## Remove duplicate rows:
synonyms_mt <- distinct(synonyms_mt)

## All results ####
all.resmt <- rbind(sp_iucnmt, syn_res_dfmt)
all.resmt <- all.resmt %>% 
  dplyr::select(taxonid, scientific_name, kingdom, phylum, class, order, family,
                genus, category, population_trend, marine_system, freshwater_system,
                terrestrial_system) %>% # keep only certain columns
  full_join(., synonyms, by = 'scientific_name') # add synonym column to match

## Merge with mtdna data ####
mtdna_iucn <- merge(mtdna, all.resmt, by.x = "SP", by.y = "synonym", all = TRUE, 
                    incomparables = NA)

which(duplicated(mtdna_iucn$SP))
View(mtdna[c(114:117, 202, 203, 205, 206, 362:365, 860:863, 931, 932, 961:963),])

# Automolus rubiginosus is Clibanornis rufipectus; keep only row 114 
# Remove Catharacta lonnbergi, subspecies of Catharacta antarctica
# Remove Falco pelegrinoides, this is a subspecies only occurring on Canary Islands
mtdna_iucn2 <- mtdna_iucn[-c(931, 362, 363, 205, 206, 115:117),]
mtdna_iucn2 <- distinct(mtdna_iucn2) # remove duplicates

## Check IUCN cateoory matches
unique(mtdna_iucn2$IUCN)
unique(mtdna_iucn2$category) ## EW = "extinct in the wild"

nomatch <- mtdna_iucn2[which(mtdna_iucn2$IUCN != mtdna_iucn2$category),]

View(nomatch %>% select(SP, IUCN, category))

mtdna_iucn2 <- mtdna_iucn2 %>% 
               mutate(IUCN_fact = as.factor(case_when(category == 'CR' ~ '5',
                                                      category == 'EN' ~ '4',
                                                      category == 'VU' ~ '3',
                                                      category == 'NT' ~ '2',
                                                      category == 'LC' ~ '1')),
               IUCN_bin = ifelse(category %in% c('CR', 'EN', 'VU'), 1, 0))
## file ------
#write.csv(mtdna_iucn2, 'mtdna_data.csv', row.names = F, quote = F)

# whole genome data ####
## Data from: Brüniche-Olsen et al. 2021 Proceedings B, Life-history traits and habitat availability shape genomic diversity in birds: implications for conservation
## https://royalsocietypublishing.org/doi/10.1098/rspb.2021.1441
wgs <- read.csv('Bruniche_Olsen_WGS_iucn.csv', h=T, na.strings = "NA")

wgs$row_id <- c(1:nrow(wgs))
wgs$to_IUCN <- paste0(wgs$genus, ' ', wgs$species)

## Search ####
species_list <- as.list(wgs$to_IUCN)

# Search over species list
sp_iucn_list <- lapply(species_list, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Create dataframe
sp_iucn <- do.call(rbind, sp_iucn_list)

## Find species with no search results ####
no_res <- which(unlist(lapply(sp_iucn_list, length))==0)
sp_not_found <- species_list[no_res]

## Find synonyms and accepted names if possible ####
syn_list <- lapply(sp_not_found, function(x) rl_synonyms(x, key = IUCN_REDLIST_KEY)$result)
syn_df <- do.call(rbind, syn_list)

# See if any had no synonyms:
no_res_syn <- which(unlist(lapply(syn_list, length))==0)
syn_not_found <- unlist(sp_not_found[no_res_syn]) 

# Create new list of species to search without duplicates
syn_search <- unique(syn_df$accepted_name)

# Search:
syn_search_res <- lapply(syn_search, function(x) rl_search(x, key = IUCN_REDLIST_KEY)$result)

# Dataframe with synonym results:
syn_res_df <- do.call(rbind, syn_search_res)

## unified taxonomy ####
synonyms <- syn_df %>% ## (!) note syn_df may have duplicate accepted names with multiple synonyms
  select(accepted_name, synonym) # synonym is the column that matches with original data

# create a dataframe with accepted & synonym to merge all results:
syn_merge <- data.frame(accepted_name = sp_iucn$scientific_name,
                        synonym = sp_iucn$scientific_name)

synonyms <- rbind(synonyms, syn_merge)
names(synonyms)[1] <- "scientific_name"

## Remove duplicate rows:
synonyms <- distinct(synonyms)

## All results ####
all.res <- rbind(sp_iucn, syn_res_df)
all.res <- all.res %>% 
  select(taxonid, scientific_name, kingdom, phylum, class, order, family,
         genus, category, population_trend, marine_system, freshwater_system,
         terrestrial_system) %>% # keep only certain columns
  full_join(., synonyms, by = 'scientific_name')

## Merge with SNP data ####
wgs_iucn <- merge(wgs, all.res, by.x = "to_IUCN", by.y = "synonym", all = TRUE, 
                  incomparables = NA)

## GD summary ####
wgs_iucn2 <- wgs_iucn %>%
  mutate(IUCN_fact = as.factor(case_when(category == 'CR' ~ '5',
                                         category == 'EN' ~ '4',
                                         category == 'VU' ~ '3',
                                         category == 'NT' ~ '2',
                                         category == 'LC' ~ '1')),
         IUCN_binary = ifelse(category %in% c('CR', 'EN', 'VU'), 
                              1, 0))

## file ####
#write.csv(wgs_iucn2, "WGS_data.csv", row.names = F, quote = F)