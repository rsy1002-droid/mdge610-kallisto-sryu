.libPaths("C:/Rlibs")
.libPaths()

setwd("D:/mdge610")
getwd()

library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(pals)
library(stringr)
library(writexl)
library(apeglm)
library(readr)


ground_truth <- read_delim("sim_true_counts.txt")

#CORRELATION (PEARSON)
##default
default <- read_tsv("abundance_dir.tsv")
default$target_id <- gsub("\\|.*", "", default$target_id)
default <- default[,c(1,4)]
colnames(default)[1] <- "transcript_id"
default <- left_join(ground_truth, default, join_by(transcript_id))
default <- default[default$true_counts >0, ]
cor.test(default$true_counts, default$est_counts, method = "pearson") # cor = 0.9925949, df = 68589, p-value < 2.2e-16
##iter1
iter1 <- read_tsv("abundance_iter_1.tsv")
iter1$target_id <- gsub("\\|.*", "", iter1$target_id)
iter1 <- iter1[,c(1,4)]
colnames(iter1)[1] <- "transcript_id"
iter1 <- left_join(ground_truth, iter1, join_by(transcript_id))
iter1 <- iter1[iter1$true_counts >0, ]
cor.test(iter1$true_counts, iter1$est_counts, method = "pearson") # cor = 0.8263167
##iter10
iter10$target_id <- gsub("\\|.*", "", iter10$target_id)
iter10 <- iter10[,c(1,4)]
colnames(iter10)[1] <- "transcript_id"
iter10 <- left_join(ground_truth, iter10, join_by(transcript_id))
iter10 <- iter10[iter10$true_counts >0, ]
cor.test(iter10$true_counts, iter10$est_counts, method = "pearson") # cor = 0.9675642
##iter25
iter25 <- read_tsv("abundance_iter_25.tsv")
iter25$target_id <- gsub("\\|.*", "", iter25$target_id)
iter25 <- iter25[,c(1,4)]
colnames(iter25)[1] <- "transcript_id"
iter25 <- left_join(ground_truth, iter25, join_by(transcript_id))
iter25 <- iter25[iter25$true_counts >0, ]
cor.test(iter25$true_counts, iter25$est_counts, method = "pearson") # cor = 0.9853814
##iter50
iter50 <- read_tsv("abundance_iter_50.tsv")
iter50$target_id <- gsub("\\|.*", "", iter50$target_id)
iter50 <- iter50[,c(1,4)]
colnames(iter50)[1] <- "transcript_id"
iter50 <- left_join(ground_truth, iter50, join_by(transcript_id))
iter50 <- iter50[iter50$true_counts >0, ]
cor.test(iter50$true_counts, iter50$est_counts, method = "pearson") # cor = 0.9907121
##iter100
iter100 <- read_tsv("abundance_iter_100.tsv")
iter100$target_id <- gsub("\\|.*", "", iter100$target_id)
iter100 <- iter100[,c(1,4)]
colnames(iter100)[1] <- "transcript_id"
iter100 <- left_join(ground_truth, iter100, join_by(transcript_id))
iter100 <- iter100[iter100$true_counts >0, ]
cor.test(iter100$true_counts, iter100$est_counts, method = "pearson") # cor = 0.9925347
##iter500
iter500 <- read_tsv("abundance_iter_500.tsv")
iter500$target_id <- gsub("\\|.*", "", iter500$target_id)
iter500 <- iter500[,c(1,4)]
colnames(iter500)[1] <- "transcript_id"
iter500 <- left_join(ground_truth, iter500, join_by(transcript_id))
iter500 <- iter500[iter500$true_counts >0, ]
cor.test(iter500$true_counts, iter500$est_counts, method = "pearson") # cor = 0.9927124
##iter1000
iter1000 <- read_tsv("abundance_iter_1000.tsv")
iter1000$target_id <- gsub("\\|.*", "", iter1000$target_id)
iter1000 <- iter1000[,c(1,4)]
colnames(iter1000)[1] <- "transcript_id"
iter1000 <- left_join(ground_truth, iter1000, join_by(transcript_id))
iter1000 <- iter1000[iter1000$true_counts >0, ]
cor.test(iter1000$true_counts, iter1000$est_counts, method = "pearson") # cor = 0.9926092
##iter2000
iter2000 <- read_tsv("abundance_iter_2000.tsv")
iter2000$target_id <- gsub("\\|.*", "", iter2000$target_id)
iter2000 <- iter2000[,c(1,4)]
colnames(iter2000)[1] <- "transcript_id"
iter2000 <- left_join(ground_truth, iter2000, join_by(transcript_id))
iter2000 <- iter2000[iter2000$true_counts >0, ]
cor.test(iter2000$true_counts, iter2000$est_counts, method = "pearson") # cor = 0.9925265
##conv.thre. 1e-1
conv1 <- read_tsv("abundance_conv_1e-1.tsv")
conv1$target_id <- gsub("\\|.*", "", conv1$target_id)
conv1 <- conv1[,c(1,4)]
colnames(conv1)[1] <- "transcript_id"
conv1 <- left_join(ground_truth, conv1, join_by(transcript_id))
conv1 <- conv1[conv1$true_counts >0, ]
cor.test(conv1$true_counts, conv1$est_counts, method = "pearson") #cor = 0.9927175
##conv.thre. 1e-3
conv3 <- read_tsv("abundance_conv_1e-3.tsv")
conv3$target_id <- gsub("\\|.*", "", conv3$target_id)
conv3 <- conv3[,c(1,4)]
colnames(conv3)[1] <- "transcript_id"
conv3 <- left_join(ground_truth, conv3, join_by(transcript_id))
conv3 <- conv3[conv3$true_counts >0, ]
cor.test(conv3$true_counts, conv3$est_counts, method = "pearson") #cor = 0.992491
##conv.thre. 1e-5
conv5 <- read_tsv("abundance_conv_1e-5.tsv")
conv5$target_id <- gsub("\\|.*", "", conv5$target_id)
conv5 <- conv5[,c(1,4)]
colnames(conv5)[1] <- "transcript_id"
conv5 <- left_join(ground_truth, conv5, join_by(transcript_id))
conv5 <- conv5[conv5$true_counts >0, ]
cor.test(conv5$true_counts, conv5$est_counts, method = "pearson") #cor = 0.9924877


#CORRELATION (PEARSON, log2 scale)
##default
cor.test(log2(default$true_counts+1), log2(default$est_counts+1), method = "pearson") # cor = 0.8007351, df = 68589, p-value < 2.2e-16
##iter1
cor.test(log2(iter1$true_counts+1), log2(iter1$est_counts+1), method = "pearson") # cor = 0.6958851 
##iter10
cor.test(log2(iter10$true_counts+1), log2(iter10$est_counts+1), method = "pearson") # cor = 0.8414484
##iter25
cor.test(log2(iter25$true_counts+1), log2(iter25$est_counts+1), method = "pearson") # cor = 0.8436043
##iter50
cor.test(log2(iter50$true_counts+1), log2(iter50$est_counts+1), method = "pearson") # cor = 0.8340164
##iter100
cor.test(log2(iter100$true_counts+1), log2(iter100$est_counts+1), method = "pearson") # cor = 0.8217908
##iter500
cor.test(log2(iter500$true_counts+1), log2(iter500$est_counts+1), method = "pearson") # cor = 0.803826
##iter1000
cor.test(log2(iter1000$true_counts+1), log2(iter1000$est_counts+1), method = "pearson") # cor = 0.8009799
##iter2000
cor.test(log2(iter2000$true_counts+1), log2(iter2000$est_counts+1), method = "pearson") # cor = 0.7994832,
##conv.thre 1e-1
cor.test(log2(conv1$true_counts+1), log2(conv1$est_counts+1), method = "pearson") # cor = 0.8187789 
##conv.thre 1e-3
cor.test(log2(conv3$true_counts+1), log2(conv3$est_counts+1), method = "pearson") # cor = 0.7983036
##conv.thre 1e-5
cor.test(log2(conv5$true_counts+1), log2(conv5$est_counts+1), method = "pearson") # cor = 0.7982354


#RELATIVE ERROR
##default
default <- default %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(default$relative_error) # Median 0.53317, Mean = 0.85207
#iter1
iter1 <- iter1 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter1$relative_error)  # Median = 0.4850, Mean = 1.8836
#iter10
iter10 <- iter10 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter10$relative_error)  # Median = 0.3779, Mean = 0.9657
#iter25
iter25 <- iter25 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter25$relative_error)  # Median = 0.40611, Mean = 0.83987
#iter50
iter50 <- iter50 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter50$relative_error) # Median = 0.44407, Mean = 0.81737
#iter100
iter100 <- iter100 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter100$relative_error) #Median = 0.48463, Mean = 0.82512
#iter500
iter500 <- iter500 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter500$relative_error) #Median = 0.5262, Mean = 0.8485
#iter1000
iter1000 <- iter1000 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter1000$relative_error) #Median = 0.53252, Mean = 0.85179
#iter2000
iter2000 <- iter2000 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(iter2000$relative_error) #Median = 0.53634, Mean = 0.85374
#conv1
conv1 <- conv1 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(conv1$relative_error) #Median = 0.49220, Mean = 0.82837
#conv3
conv3 <- conv3 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(conv3$relative_error) #Median = 0.53832, Mean = 0.85508
#conv5
conv5 <- conv5 %>% mutate(relative_error = abs(est_counts - true_counts) / true_counts)
summary(conv5$relative_error) #Median = 0.53838, Mean = 0.85511


#make a dataframe
accuracy <- data.frame(
  setting = c("default", "iter 1", "iter 10", "iter 25", "iter 50", "iter 100", "iter 500", "iter 1000", "iter 2000", "thre 1e-1", "thre 1e-3", "thre 1e-5"),
  iteration_number = c(1095, 1, 10, 25, 50, 100, 500, 1000, 2000, 119, 7828, 10000),
  pearson_log2 = c(0.8007351, 0.6958851, 0.8414484, 0.8436043, 0.8340164, 0.8217908, 0.803826, 0.8009799, 0.7994832, 0.8187789, 0.7983036, 0.7982354),
  relative_error_median = c(0.53317, 0.4850, 0.3779, 0.40611, 0.44407, 0.48463, 0.5262, 0.53252, 0.53634, 0.49220, 0.53832, 0.53838),
  relative_error_mean = c(0.85207, 1.8836, 0.9657, 0.83987, 0.81737, 0.82512, 0.8485, 0.85179, 0.85374, 0.82837, 0.85508, 0.85511)
)


#make plots 
accuracy <- accuracy[order(accuracy$iteration_number), ]

p1 <- ggplot(accuracy, aes(x = iteration_number, y = pearson_log2))+
  geom_line()+
  scale_x_log10()+
  labs(title = "Accuracy Analysis: Pearson Correlation vs. Iteration Number", x = "iteration_number (log scale)", y = "pearson_log2")+
  theme_minimal()

p2 <- ggplot(accuracy, aes(x = iteration_number, y = relative_error_median))+
  geom_line()+
  scale_x_log10()+
  labs(title = "Accuracy Analysis: Median Relative Error vs. Iteration Number", x = "iteration_number (log scale)", y = "relative_error_median")+
  theme_minimal()


ggarrange(p1, p2, ncol = 2)

