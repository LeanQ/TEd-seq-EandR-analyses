library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(patchwork)
library(gridExtra)
library(reshape2)
library("ggsci")

# We simulate each generation with the following code;


simulate_G1<-function(p) {
N=2500 #population size at generation 1
fAA_in <- p^2
fAa_in <- 2*p*(1-p)
faa_in <- (1-p)^2

s <- 0.95 #selfing rate

fAA <- s*(fAA_in + fAa_in/4) + (1-s)*p^2
fAa <- s*fAa_in/2 + (1-s)*2*p*(1-p)
faa <- s*(faa_in + fAa_in/4) + (1-s)*(1-p)^2
p <- fAA + faa/2

fAA
fAa
AA <- rbinom(1, N, fAA)
Aa <- rbinom(1, N, fAa)
AA
Aa
A_alleles<-2*AA+Aa
total_alleles <-2*N
p_next <-A_alleles/total_alleles
return(p_next)
p_next
}

simulate_G2<-function(p) {
N=9000 #population size 
fAA_in <- p^2
fAa_in <- 2*p*(1-p)
faa_in <- (1-p)^2

s <- 0.95 #selfing rate

fAA <- s*(fAA_in + fAa_in/4) + (1-s)*p^2
fAa <- s*fAa_in/2 + (1-s)*2*p*(1-p)
faa <- s*(faa_in + fAa_in/4) + (1-s)*(1-p)^2
p <- fAA + faa/2

AA <- rbinom(1, N, fAA)
Aa <- rbinom(1, N, fAa)
A_alleles<-2*AA+Aa
total_alleles <-2*N
p_next<-A_alleles/total_alleles
return (p_next)

}

set.seed(32) 
# We estimate the mean frequency and standard deviation taking into acount the replicates
calculate_statistics <- function(observed_freq) {
  mean_observed <- mean(observed_freq)
  sd_observed <- sd(observed_freq)
  return(list(mean_observed = mean_observed, sd_observed = sd_observed))
}


#analysis for population 1

population_1 <- read.table("/path_to_file/allele_freq_population1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

unique_G0 <- data.frame() #insertions that are specific of G0, will be discarded. Probably specific of the pool of 1000 seeds 
unique_G1 <- data.frame() #insertions that are specific of G1, will be discarded. Probably specific of the pools of 1000 seeds 
unique_G2 <- data.frame() #insertions that are specific of G2, will be discarded. Probably specific of the pools of 1000 seeds 
new_insertions_G1 <- data.frame() # to store all those insertions that are novel of G1 and detected in generation 1
filtered_shared_insertions <- data.frame() #insertions heritable that we can detect across generations and have not appear in generation 1, those are the ones that we will proceed to analyze

# Function to classify rows
classify_rows <- function(row) {
  # Convert relevant columns to numeric
  G0R1 <- as.numeric(row[['G0R1']])
  G0R2 <- as.numeric(row[['G0R2']])
  G1P1R1 <- as.numeric(row[['G1P1R1']])
  G1P1R2 <- as.numeric(row[['G1P1R2']])
  G1P1R3 <- as.numeric(row[['G1P1R3']])
  G2P1R1 <- as.numeric(row[['G2P1R1']])
  G2P1R2 <- as.numeric(row[['G2P1R2']])
  G2P1R3 <- as.numeric(row[['G2P1R3']])
  
  # Check conditions and classify
  if (all(c(G0R1, G0R2, G1P1R1, G1P1R2, G1P1R3, G2P1R1, G2P1R2, G2P1R3) < 0.001)) {
    return("discard")
  }
  
  if (any(c(G0R1, G0R2) > 0.001) && all(c(G1P1R1, G1P1R2, G1P1R3, G2P1R1, G2P1R2, G2P1R3) < 0.001)) {
    return("unique_G0")
  }
  
  if (any(c(G1P1R1, G1P1R2, G1P1R3) > 0.001) && all(c(G0R1, G0R2, G2P1R1, G2P1R2, G2P1R3) < 0.001)) {
    return("unique_G1")
  }
  
  if (any(c(G2P1R1, G2P1R2, G2P1R3) > 0.001) && all(c(G0R1, G0R2, G1P1R1, G1P1R2, G1P1R3) < 0.001)) {
    return("unique_G2")
  }
  
  if (any(c(G1P1R1, G1P1R2, G1P1R3) > 0.001) && any(c(G2P1R1, G2P1R2, G2P1R3) > 0.001) && all(c(G0R1, G0R2) < 0.001)) {
    return("new_insertions_G1")
  }
  
  return("filtered_shared_insertions")
}

# Apply classification function to each row and split into different dataframes
population_1$classification <- apply(population_1, 1, classify_rows)

pop1_unique_G0 <- population_1[population_1$classification == "unique_G0", ]
pop1_unique_G1 <- population_1[population_1$classification == "unique_G1", ]
pop1_unique_G2 <- population_1[population_1$classification == "unique_G2", ]
pop1_new_insertions_G1 <- population_1[population_1$classification == "new_insertions_G1", ]
pop1_shared_insertions <- population_1[population_1$classification == "filtered_shared_insertions", ]


#We proceed to filter out all those insertions that are at a very low coverage

frequency_threshold <- 0.002 #all insertions lower that have a lower frequency will not be considered
filtered_shared_insertions_1 <- pop1_shared_insertions[
  apply(pop1_shared_insertions[, c('G0R1', 'G0R2', 'G1P1R1', 'G1P1R2', 'G1P1R3', 'G2P1R1', 'G2P1R2', 'G2P1R3')], 1, 
        function(row) any(row > frequency_threshold, na.rm = TRUE)), 
]

# Then we proceed to obtain the observed allele frequencies and compare to the expected by genetic drift
for (i in 1:nrow(filtered_shared_insertions_1)) {
  row <- filtered_shared_insertions_1[i, ]
  
  # We extract observed frequencies
  observed_freq_gen0 <- c(row[['G0R1']], row[['G0R2']])
  observed_freq_gen1 <- c(row[['G1P1R1']], row[['G1P1R2']], row[['G1P1R3']])
  observed_freq_gen2 <- c(row[['G2P1R1']], row[['G2P1R2']], row[['G2P1R3']])
  
  #We calculate statistics for observed frequencies
  stats_gen0 <- calculate_statistics(observed_freq_gen0)
  stats_gen1 <- calculate_statistics(observed_freq_gen1)
  stats_gen2 <- calculate_statistics(observed_freq_gen2)
  
  # Initial frequency estimates (G0R1 and G0R2)
  initial_frequencies <- c(as.numeric(row[['G0R1']]), as.numeric(row[['G0R2']]))
  
  # Calculate mean initial frequency
  p_initial <- mean(initial_frequencies)
  
  #We simulate expected allelic frequency for generation 1 and generation 2. We do for each insertion 1000 simulations
  
  expected_gen1 <- replicate(1000, simulate_G1( p_initial))
  expected_gen1_mean <- mean(expected_gen1)
  expected_gen1_sd <- sd(expected_gen1)
  
  expected_gen2 <- replicate(1000, simulate_G2( expected_gen1_mean))
  expected_gen2_mean <- mean(expected_gen2)
  expected_gen2_sd <- sd(expected_gen2)
  
  #We obtain a p-value for observed vs. expected allelic frequency in generation 2 by comparing using a zscore of the observed vs expected allele frequencies
  n <- length(observed_freq_gen2)
  
  z_score <- (stats_gen2$mean_observed - expected_gen2_mean) / (stats_gen2$sd_observed / sqrt(n))
    p_value <- 2 * pnorm(-abs(z_score))  # Two-tailed test
  # We classify finally the insertions based on the pattern observed and the pvalue
    
    if (stats_gen1$mean_observed > stats_gen0$mean_observed && stats_gen2$mean_observed > stats_gen1$mean_observed) {
      classification <- "Increased"
    } else if ( stats_gen2$mean_observed < stats_gen1$mean_observed &&stats_gen1$mean_observed < stats_gen0$mean_observed ) {
      classification <- "Decreased"
    } else {
      classification <- "Not_sig"
    }
      filtered_shared_insertions_1[i, 'p_value'] <- p_value
      filtered_shared_insertions_1[i, 'classification'] <- classification

  } 
 #Finally we correct the pvalues by the false discovery test and adjust the data accordingly
filtered_shared_insertions_1['padj']<- p.adjust(filtered_shared_insertions_1$p_value,method="fdr")

filtered_shared_insertions_1<-filtered_shared_insertions_1 %>%
  mutate(classification = if_else(padj > 0.05, "Not_sig", classification))
View(filtered_shared_inserti)
for (i in 1:nrow(filtered_shared_insertions_1)) {
  row <- filtered_shared_insertions_1[i, ]
 if (row$padj >0.05 ) {
      classification <- "Not_sig"
}}
