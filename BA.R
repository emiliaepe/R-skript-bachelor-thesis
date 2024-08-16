################################################################################
############################### Bachelor thesis ################################
################# Characterization of mucin-degrading bacteria #################
################################################################################

# data import
install.packages("readr")
library(readr)

BA_R <- read_delim("Desktop/BA_R_neu.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(BA_R)

# separate by database
data_ltp <- subset(BA_R, BA_R$data_base == "LTP")
data_ncbi <- subset(BA_R, BA_R$data_base == "NCBI")
str(data_ncbi)

######################## 3.2 Comparison of the sites ##########################
######################## Composition of phyla plot ############################ 

# determine the order of the locations
data_ncbi$location <- factor(data_ncbi$location, levels = c("Svalbard", "North Sea", "Molgula"))

# phyla per location
install.packages("dplyr")
library(dplyr)

phyla_counts <- data_ncbi %>%
  group_by(location, phylum) %>%
  summarise(count = n()) %>%
  ungroup()

phyla_totals <- phyla_counts %>%
  group_by(location) %>%
  summarise(total = sum(count))

phyla_percent <- phyla_counts %>%
  left_join(phyla_totals, by = "location") %>%
  mutate(percent = (count / total) * 100)

# create a stacked bar chart
install.packages("ggplot2")
install.packages("scico")
library(ggplot2)
library(scico)

ggplot(phyla_percent, aes(x = location, y = percent, fill = phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_scico_d(palette = "roma") +
  labs(title = "Site-specific composition of Phyla",
       x = "Site",
       y = "Percentage composition",
       fill = "Phylum") +
  theme_minimal()

############################ additional plots ################################# 

# Ffunction for stacked bar chart
create_stacked_bar_plot <- function(data, tax_level) {
  tax_counts <- data %>%
    group_by(location, .data[[tax_level]]) %>%
    summarise(count = n()) %>%
    ungroup()
  
  tax_totals <- tax_counts %>%
    group_by(location) %>%
    summarise(total = sum(count))
  
  tax_percent <- tax_counts %>%
    left_join(tax_totals, by = "location") %>%
    mutate(percent = (count / total) * 100)
  
  ggplot(tax_percent, aes(x = location, y = percent, fill = .data[[tax_level]])) +
    geom_bar(stat = "identity") +
    scale_fill_scico_d(palette = "roma") +
    labs(title = paste("Prozentuale Zusammensetzung der", tax_level, "pro Standort"),
         x = "Standort",
         y = "Prozentuale Zusammensetzung",
         fill = tax_level) +
    theme_minimal()
}

# create diagrams for various tax. level 
plot_phylum <- create_stacked_bar_plot(data_ncbi, "phylum")
plot_family <- create_stacked_bar_plot(data_ncbi, "family")
plot_genus <- create_stacked_bar_plot(data_ncbi, "genus")

# show diagram
print(plot_phylum)
print(plot_family)
print(plot_genus)

########################### site-specific plots ###############################
install.packages("plotly")
library(plotly)

# separate by location
data_north_sea <- data_ncbi %>%
  filter(location == "North Sea") 

data_molgula <-  data_ncbi %>%
  filter(location == "Molgula") 

data_svalbard <-  data_ncbi %>%
  filter(location == "Svalbard")
 
# define color for family
unique_families <- unique(data_ncbi$family)
family_colors <- scico::scico(length(unique_families), palette = "roma")  # Wähle eine Palette wie "berlin" oder "lajolla"
names(family_colors) <- unique_families

# function for creating the nested pie chart
plot_nested_pie_chart <- function(data, dataset_name) {
  
  # aggregate data by family and genus
  genus_data <- data %>%
    group_by(family, genus) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(percentage = count / sum(count) * 100)
  
  # aggregate for families
  family_data <- genus_data %>%
    group_by(family) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    mutate(percentage = count / sum(count) * 100)
  
  # order in family_colors
  family_data$family <- factor(family_data$family, levels = names(family_colors))
  genus_data$family <- factor(genus_data$family, levels = names(family_colors))
  
  # sort by the factor levels
  family_data <- family_data[order(family_data$family), ]
  genus_data <- genus_data[order(genus_data$family), ]
  
  # Nested pie chart
  plot_ly() %>%
    add_pie(
      data = family_data,
      labels = ~family,
      values = ~count,
      hole = 0.6,
      name = "Family",
      textinfo = 'label+percent',
      sort = FALSE,  # not sort by value
      marker = list(colors = family_colors[family_data$family]),
      legendgroup = 'family'
    ) %>%
    add_pie(
      data = genus_data,
      labels = ~genus,
      values = ~count,
      hole = 0.4,
      name = "Genus",
      textinfo = 'label+percent',
      sort = FALSE,  # not sort by value
      marker = list(
        colors = family_colors[genus_data$family]
      ),
      showlegend = FALSE
    ) %>%
    layout(
     # title = paste("Nested Pie Chart der Bakteriengemeinschaft für", dataset_name),
      showlegend = TRUE,
      legend = list(
        x = 0.8,
        y = 0.5,
        traceorder = 'normal',
        orientation = 'v'
      )
    )
}

# show diagram
plot_nested_pie_chart(data_svalbard, "Svalbard")
plot_nested_pie_chart(data_molgula, "Molgula sp.")
plot_nested_pie_chart(data_north_sea, "North Sea")

########################## Abundancy calculations #############################

# calculation of the frequencies of genera and families
genus_counts <- data_molgula %>%
  count(genus, family) %>%
  rename(count = n)

# calculation of the total number of genera and families
total_count <- sum(genus_counts$count)

# calculation of the percentage of genera
genus_counts <- genus_counts %>%
  mutate(percent_genus = round(100 * count / total_count, 1))

# frequency table for genus
taxon_matrix_ltp <- table(data_ltp$genus)
taxon_matrix_ncbi <- table(data_ncbi$genus)

taxon_matrix_north_sea <- table(data_north_sea$genus)
taxon_matrix_molgula <- table(data_molgula$genus)
taxon_matrix_svalbard <- table(data_svalbard$genus)

# create contingency table
contingency_table <- table(data_ncbi$location, data_ncbi$genus)

# show table
print(contingency_table)
View(contingency_table)
str(contingency_table)

# convert contingency table to DataFrame
contingency_df <- as.data.frame(contingency_table)

# define order of locations
contingency_df$Var1 <- factor(contingency_df$Var1, levels = c("Svalbard","North Sea", "Molgula"))

# calculate number of different genera per location
genus_count_per_location <- data_ncbi %>%
  group_by(location) %>%
  summarize(num_genus = n_distinct(genus))

# show result
print(genus_count_per_location)

# only in 1 Location
single_location_genus <- contingency_df %>%
  group_by(Var2) %>%
  filter(sum(Freq > 0) == 1) %>%  # Genus muss in genau einer Location vorkommen
  filter(Freq > 0) %>%             # Stelle sicher, dass die Frequenz in dieser Location > 0 ist
  select(Var2, Var1) %>%
  ungroup()

# show result 
View(single_location_genus)

# export as CSV file
write.csv(single_location_genus, file = "single_location_genus.csv", row.names = FALSE)

# in 2 location 
two_location_genus <- contingency_df %>%
  group_by(Var2) %>%
  filter(sum(Freq > 0) == 2) %>%  # Genus muss in genau einer Location vorkommen
  filter(Freq > 0) %>%             # Stelle sicher, dass die Frequenz in dieser Location > 0 ist
  select(Var2, Var1) %>%
  ungroup()

# show result 
View(two_location_genus)

# export as CSV file
write.csv(two_location_genus, file = "two_location_genus.csv", row.names = FALSE)

# in 3 locations 
three_location_genus <- contingency_df %>%
  group_by(Var2) %>%
  filter(sum(Freq > 0) == 3) %>%  # Genus muss in genau einer Location vorkommen
  filter(Freq > 0) %>%             # Stelle sicher, dass die Frequenz in dieser Location > 0 ist
  select(Var2, Var1) %>%
  ungroup()

# show result 
View(three_location_genus)

# export as CSV file
write.csv(three_location_genus, file = "three_location_genus.csv", row.names = FALSE)

########################### Heatap Genus Abundance ############################

# create Heatmap
ggplot(contingency_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#00BFC4") +
  labs(x = "Site", y = "Genus", fill = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


########################### alpha-diversitys sites #############################

shannon_index_north_sea <- diversity(taxon_matrix_north_sea, index = "shannon")
simpson_index_north_sea <- diversity(taxon_matrix_north_sea, index = "simpson")

shannon_index_molgula <- diversity(taxon_matrix_molgula, index = "shannon")
simpson_index_molgula <- diversity(taxon_matrix_molgula, index = "simpson")

shannon_index_svalbard <- diversity(taxon_matrix_svalbard, index = "shannon")
simpson_index_svalbard <- diversity(taxon_matrix_svalbard, index = "simpson")


# data frame for the diversity indices
diversity_summary_location <- data.frame(
  Site = c("North Sea", "Molgula sp.", "Svalbard"),
  Shannon_Index = c(shannon_index_north_sea, shannon_index_molgula, shannon_index_svalbard),
  Simpson_Index = c(simpson_index_north_sea, simpson_index_molgula, simpson_index_svalbard)
)

print(diversity_summary_location)

#  visualization
diversity_df_location <- data.frame(
  Group = rep(c("Svalbard", "North Sea", "Molgula sp."), each = 2),
  Index = c(shannon_index_svalbard, simpson_index_svalbard, shannon_index_north_sea, 
            simpson_index_north_sea, shannon_index_molgula, simpson_index_molgula),
  Type = rep(c("Shannon", "Simpson"), times = 3)
)

# set the order of the groups
diversity_df_location$Group <- factor(diversity_df_location$Group, levels = c("Svalbard", "North Sea", "Molgula sp."))

# create the plot
p <- ggplot(diversity_df_location, aes(x = Group, y = Index, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#7E1700", "#00BFC4", "#00BFC4")) +
  theme_minimal() +
  labs(title = "Comparison of the sites α diversity indices", x = "Site", y = "Index value")

ggplot_build(p)$data[[1]]$fill %>% unique()
"#7E1700" "#00BFC4"

########################### beta-diversity sites #############################
install.packages("tidyr")
library(tidyr)


str(contingency_df)
contingency_wide <- spread(contingency_df, key = Var2, value = Freq, fill = 0)
View(contingency_wide)
contingency_wide <- pivot_wider(contingency_df, names_from = Var2, values_from = Freq, values_fill = list(Freq = 0))


# converting the contingency table into a numerical matrix
contingency_matrix <- as.matrix(contingency_wide[,-1]) # Entfernen der ersten Spalte, die die Standort-Namen enthält
rownames(contingency_matrix) <- as.character(contingency_wide$Var1) # Setzen der Standortnamen als Reihenbezeichner
View(contingency_matrix)

# convert the matrix values into numerical values
contingency_matrix <- apply(contingency_matrix, 2, as.numeric)
View(contingency_matrix)
str(contingency_matrix)

# calculating the Bray-Curtis dissimilarity
bray_curtis <- vegdist(contingency_matrix, method = "bray")

# output of the Bray-Curtis dissimilarity
print(bray_curtis)

# output in data frame
bray_curtis_df <- as.data.frame(as.matrix(bray_curtis))
View(bray_curtis_df)


#################### 4.4.4	Comparison of the databases  ######################
######################## alpha-diversity databases ############################

# contingency table data comparison
contingency_table <- read_delim("Desktop/Kontingenztabelle_genera_ databases.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
# exact Fisher Test 
fisher_test <- fisher.test(contingency_table[, 2:3])
print(fisher_test)

# calculation of Shannon and Simpson index
install.packages("vegan")
library(vegan)

shannon_index_ltp <- diversity(taxon_matrix_ltp, index = "shannon")
simpson_index_ltp <- diversity(taxon_matrix_ltp, index = "simpson")

shannon_index_ncbi <- diversity(taxon_matrix_ncbi, index = "shannon")
simpson_index_ncbi <- diversity(taxon_matrix_ncbi, index = "simpson")

# data frame for indices 
diversity_summary_database <- data.frame(
  Database = c("LTP", "NCBI"),
  Shannon_Index = c(shannon_index_ltp, shannon_index_ncbi),
  Simpson_Index = c(simpson_index_ltp, simpson_index_ncbi)
)

print(diversity_summary_database)

#  visualization 
diversity_df_database <- data.frame(
  Group = rep(c("LTP", "NCBI"), each = 2),
  Index = c(shannon_index_ltp, simpson_index_ltp,
            shannon_index_ncbi, simpson_index_ncbi),
  Type = rep(c("Shannon", "Simpson"), times = 2)
)

# create a bar chart
ggplot(diversity_df_database, aes(x = Group, y = Index, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#7E1700", "#00BFC4", "#00BFC4"))+
  theme_minimal() +
  labs(title = "Comparison of the α-diversity indices of the databases at genus level", x = "Database", y = "Index value")
  
######################### beta-Diversity databases ############################

# find unique genera in each data frame
genus_ltp <- unique(data_ltp$genus)
genus_ncbi <- unique(data_ncbi$genus)

# calculate intersection and union of the genera
intersection_genus <- length(intersect(genus_ltp, genus_ncbi))
union_genus <- length(union(genus_ltp, genus_ncbi))

# calculate Jaccard index
jaccard_index <- intersection_genus / union_genus

# output
print(jaccard_index)

########################### additional calculations ############################

# calculate genera only present in one database
only_in_ltp <- setdiff(genus_ltp, genus_ncbi)
only_in_ncbi <- setdiff(genus_ncbi, genus_ltp)

# number of genera that are only present in one database
only_in_ltp_count <- length(only_in_ltp)
only_in_ncbi_count <- length(only_in_ncbi)

# total number of genera in the union
total_union <- length(union_genus)

# calculate the fraction of genera that are only present in one database
fraction_only_in_ltp <- only_in_ltp_count / total_union
fraction_only_in_ncbi <- only_in_ncbi_count / total_union

# output of the results
print(paste("Gattungen nur in LTP:", only_in_ltp_count, "(", round(fraction_only_in_ltp * 100, 2), "%)"))
print(paste("Gattungen nur in NCBI:", only_in_ncbi_count, "(", round(fraction_only_in_ncbi * 100, 2), "%)"))

# output of the genera only available in one database
cat("Gattungen nur in LTP:\n")
print(only_in_ltp)

cat("\nGattungen nur in NCBI:\n")
print(only_in_ncbi)

# calculate the thirds that do not match
total_diff <- only_in_ltp_count + only_in_ncbi_count
fraction_diff <- total_diff / total_union

selected_genera <- c("Altererythrobacter", "Wocania", "Sphingorhabdus", "Lutimonas", 
                     "Rhodopirellula", "Halocynthiibacter", "Verrucomicrobium", 
                     "Planktotalea", "Tsuneonella", "Pontixanthobacter", 
                     "Cognaticolwellia", "Pseudalgibacter", "Parasphingorhabdus", 
                     "Allorhodopirellula", "Falsihalocynthiibacter", "Aliiroseovarius")

# filter rows that are contained in the list of selected genera
filtered_df <- BA_R %>% filter(genus %in% selected_genera)

# output 
print(filtered_df)
View(filtered_df)
