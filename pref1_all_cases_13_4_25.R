################ read and process csv files of the 6 cases


setwd("C:/Users/haymanay/Documents/Aya/Speciation")


# number of runs (different seeds) in each setting
sample_size <- 10

num_of_individuals_checked <- 100

# bootstrap library
library(boot)

# plot CIs library
library(plotrix)

# library for U test (Mann-Whitney-Wilcoxon)
library(coin)

# library for exact U test (Mann-Whitney-Wilcoxon) with ties
library(exactRankTests)

# library for plotting
library(ggplot2)

# library for boxplots
#library(Rlab)
library(fields)

# library for list of all possible permutations
library(permute)

# library for Levene's test for homoscedasticity
library(car)




######## tables of number of species in each checked generation for each case.
######## Each line is a different checked generation.
######## Each column is for a different run (seed).



#### table for the regular BDM case

# number of generations to look at
 num_gen_look <- 99

number_of_species_table_BDM1 <- read.table("number_of_species BDM seed_0.csv", header=FALSE)  
number_of_species_table_BDM2 <- read.table("number_of_species BDM seed_2.csv", header=FALSE)
number_of_species_table_BDM3 <- read.table("number_of_species BDM seed_3.csv", header=FALSE)  
number_of_species_table_BDM4 <- read.table("number_of_species BDM seed_4.csv", header=FALSE)
number_of_species_table_BDM5 <- read.table("number_of_species BDM seed_5.csv", header=FALSE)  
number_of_species_table_BDM6 <- read.table("number_of_species BDM seed_6.csv", header=FALSE)
number_of_species_table_BDM7 <- read.table("number_of_species BDM seed_7.csv", header=FALSE)  
number_of_species_table_BDM8 <- read.table("number_of_species BDM seed_8.csv", header=FALSE)
number_of_species_table_BDM9 <- read.table("number_of_species BDM seed_9.csv", header=FALSE)  
number_of_species_table_BDM10 <- read.table("number_of_species BDM seed_10.csv", header=FALSE)

number_of_species_table_BDM <- cbind(number_of_species_table_BDM1,
                                  number_of_species_table_BDM2,
                                  number_of_species_table_BDM3,
                                  number_of_species_table_BDM4,
                                  number_of_species_table_BDM5,
                                  number_of_species_table_BDM6,
                                  number_of_species_table_BDM7,
                                  number_of_species_table_BDM8,
                                  number_of_species_table_BDM9,
                                  number_of_species_table_BDM10)[1:num_gen_look, ]


number_of_species_table_BDM[1:5,]

(num_of_generations_checked__BDM <- dim(number_of_species_table_BDM)[1])

# vector of checked generations
checked_generations_vector_BDM <- seq(999, 98999, 1000)

# Add the fact that in the begining there is one species
checked_generation_vector_BDM_with_gen_0 <- c(0, checked_generations_vector_BDM)
number_of_species_table_BDM_with_gen_0 <- rbind(rep(1, sample_size), number_of_species_table_BDM)



checked_generations_vector_BDM_0 <- checked_generations_vector_BDM

# vector of the mean number of species for each checked generation
mean_num_of_species_vector_BDM <- apply(X=number_of_species_table_BDM, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_BDM <- apply(X=number_of_species_table_BDM, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_BDM <- apply(X=number_of_species_table_BDM, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_BDM <- apply(X=number_of_species_table_BDM, MARGIN=1, FUN=max)




# plot number of species in each run
# log scale on y axis with range
# (plot all_BDM_runs_log_range.png, plot all_BDM_runs_log_range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_BDM))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1))
y_tic_marks_values <- c(1, 2, 5)
lo <- loess((log(median_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='All types', cex.main=0.85, xlim=c(0, 0000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(median_num_of_species_vector_BDM + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold assortative mating & Haploid threshold assortative mating', cex.main=0.75, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(upper_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(lower_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
par(new=T)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_50), y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_100 <- c(0, checked_generations_vector_BDM + 100)
par(new=T)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_100), y=log(median_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_150 <- c(0, checked_generations_vector_BDM + 150)
par(new=T)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_200 <- c(0, (checked_generations_vector_BDM + 200)[1:5])
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(median_num_of_species_vector_Haploid_threshold) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(upper_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(lower_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_250 <- c(0, (checked_generations_vector_BDM + 250)[1:5])
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_150), y=log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_250), y=log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_250), y=log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
legend("topright", legend=c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating", "Haploid threshold", "Haploid threshold assortative mating"), col=c("black", "blue", "red", "purple", "orange", "green"), fill=c("black", "blue", "red", "purple", "orange", "green"))





#### table for BDM assortative mating

number_of_species_table_BDM_assortative_mating1 <- read.table("number_of_species BDM_assortative_mating seed_0.csv", header=FALSE)  
number_of_species_table_BDM_assortative_mating2 <- read.table("number_of_species BDM_assortative_mating seed_2.csv", header=FALSE)
number_of_species_table_BDM_assortative_mating3 <- read.table("number_of_species BDM_assortative_mating seed_3.csv", header=FALSE)  
number_of_species_table_BDM_assortative_mating4 <- read.table("number_of_species BDM_assortative_mating seed_4.csv", header=FALSE)
number_of_species_table_BDM_assortative_mating5 <- read.table("number_of_species BDM_assortative_mating seed_5.csv", header=FALSE)  
number_of_species_table_BDM_assortative_mating6 <- read.table("number_of_species BDM_assortative_mating seed_6.csv", header=FALSE)
number_of_species_table_BDM_assortative_mating7 <- read.table("number_of_species BDM_assortative_mating seed_7.csv", header=FALSE)  
number_of_species_table_BDM_assortative_mating8 <- read.table("number_of_species BDM_assortative_mating seed_8.csv", header=FALSE)
number_of_species_table_BDM_assortative_mating9 <- read.table("number_of_species BDM_assortative_mating seed_9.csv", header=FALSE)  
number_of_species_table_BDM_assortative_mating10 <- read.table("number_of_species BDM_assortative_mating seed_10.csv", header=FALSE)

number_of_species_table_BDM_assortative_mating <- cbind(number_of_species_table_BDM_assortative_mating1,
                                                      number_of_species_table_BDM_assortative_mating2,
                                                      number_of_species_table_BDM_assortative_mating3,
                                                      number_of_species_table_BDM_assortative_mating4,
                                                      number_of_species_table_BDM_assortative_mating5,
                                                      number_of_species_table_BDM_assortative_mating6,
                                                      number_of_species_table_BDM_assortative_mating7,
                                                      number_of_species_table_BDM_assortative_mating8,
                                                      number_of_species_table_BDM_assortative_mating9,
                                                      number_of_species_table_BDM_assortative_mating10)[1:num_gen_look, ]

number_of_species_table_BDM_assortative_mating[1:10,]

(num_of_generations_checked__BDM_assortative_mating <- dim(number_of_species_table_BDM_assortative_mating)[1])

# vector of checked generations
checked_generations_vector_BDM_assortative_mating <- seq(999, 98999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_BDM_assortative_mating <- apply(X=number_of_species_table_BDM_assortative_mating, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_BDM_assortative_mating <- apply(X=number_of_species_table_BDM_assortative_mating, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_BDM_assortative_mating <- apply(X=number_of_species_table_BDM_assortative_mating, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_BDM_assortative_mating <- apply(X=number_of_species_table_BDM_assortative_mating, MARGIN=1, FUN=max)



#### table for Diploid threshold

number_of_species_table_Diploid_threshold1 <- read.table("number_of_species threshold seed_0.csv", header=FALSE)  
number_of_species_table_Diploid_threshold2 <- read.table("number_of_species threshold seed_2.csv", header=FALSE)
number_of_species_table_Diploid_threshold3 <- read.table("number_of_species threshold seed_3.csv", header=FALSE)  
number_of_species_table_Diploid_threshold4 <- read.table("number_of_species threshold seed_4.csv", header=FALSE)
number_of_species_table_Diploid_threshold5 <- read.table("number_of_species threshold seed_5.csv", header=FALSE)  
number_of_species_table_Diploid_threshold6 <- read.table("number_of_species threshold seed_6.csv", header=FALSE)
number_of_species_table_Diploid_threshold7 <- read.table("number_of_species threshold seed_7.csv", header=FALSE)  
number_of_species_table_Diploid_threshold8 <- read.table("number_of_species threshold seed_8.csv", header=FALSE)
number_of_species_table_Diploid_threshold9 <- read.table("number_of_species threshold seed_9.csv", header=FALSE)  
number_of_species_table_Diploid_threshold10 <- read.table("number_of_species threshold seed_10.csv", header=FALSE)

number_of_species_table_Diploid_threshold <- cbind(number_of_species_table_Diploid_threshold1,
                                               number_of_species_table_Diploid_threshold2,
                                               number_of_species_table_Diploid_threshold3,
                                               number_of_species_table_Diploid_threshold4,
                                               number_of_species_table_Diploid_threshold5,
                                               number_of_species_table_Diploid_threshold6,
                                               number_of_species_table_Diploid_threshold7,
                                               number_of_species_table_Diploid_threshold8,
                                               number_of_species_table_Diploid_threshold9,
                                               number_of_species_table_Diploid_threshold10)


number_of_species_table_Diploid_threshold[1:5,]

(num_of_generations_checked__Diploid_threshold <- dim(number_of_species_table_Diploid_threshold)[1])

# vector of checked generations
checked_generations_vector_Diploid_threshold <- seq(999, 148999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_Diploid_threshold <- apply(X=number_of_species_table_Diploid_threshold, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_Diploid_threshold <- apply(X=number_of_species_table_Diploid_threshold, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_Diploid_threshold <- apply(X=number_of_species_table_Diploid_threshold, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_Diploid_threshold <- apply(X=number_of_species_table_Diploid_threshold, MARGIN=1, FUN=max)



#### table for Diploid threshold assortative mating

number_of_species_table_Diploid_threshold_assortative_mating1 <- read.table("number_of_species threshold assortative_mating seed_0.csv", header=FALSE)  
number_of_species_table_Diploid_threshold_assortative_mating2 <- read.table("number_of_species threshold assortative_mating seed_2.csv", header=FALSE)
number_of_species_table_Diploid_threshold_assortative_mating3 <- read.table("number_of_species threshold assortative_mating seed_3.csv", header=FALSE)  
number_of_species_table_Diploid_threshold_assortative_mating4 <- read.table("number_of_species threshold assortative_mating seed_4.csv", header=FALSE)
number_of_species_table_Diploid_threshold_assortative_mating5 <- read.table("number_of_species threshold assortative_mating seed_5.csv", header=FALSE)  
number_of_species_table_Diploid_threshold_assortative_mating6 <- read.table("number_of_species threshold assortative_mating seed_6.csv", header=FALSE)
number_of_species_table_Diploid_threshold_assortative_mating7 <- read.table("number_of_species threshold assortative_mating seed_7.csv", header=FALSE)  
number_of_species_table_Diploid_threshold_assortative_mating8 <- read.table("number_of_species threshold assortative_mating seed_8.csv", header=FALSE)
number_of_species_table_Diploid_threshold_assortative_mating9 <- read.table("number_of_species threshold assortative_mating seed_9.csv", header=FALSE)  
number_of_species_table_Diploid_threshold_assortative_mating10 <- read.table("number_of_species threshold assortative_mating seed_10.csv", header=FALSE)

number_of_species_table_Diploid_threshold_assortative_mating <- cbind(number_of_species_table_Diploid_threshold_assortative_mating1,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating2,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating3,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating4,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating5,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating6,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating7,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating8,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating9,
                                                                  number_of_species_table_Diploid_threshold_assortative_mating10)

number_of_species_table_Diploid_threshold_assortative_mating[1:5,]

(num_of_generations_checked__Diploid_threshold_assortative_mating <- dim(number_of_species_table_Diploid_threshold_assortative_mating)[1])

# vector of checked generations
checked_generations_vector_Diploid_threshold_assortative_mating <- seq(999, 148999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_Diploid_threshold_assortative_mating <- apply(X=number_of_species_table_Diploid_threshold_assortative_mating, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_Diploid_threshold_assortative_mating <- apply(X=number_of_species_table_Diploid_threshold_assortative_mating, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_Diploid_threshold_assortative_mating <- apply(X=number_of_species_table_Diploid_threshold_assortative_mating, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_Diploid_threshold_assortative_mating <- apply(X=number_of_species_table_Diploid_threshold_assortative_mating, MARGIN=1, FUN=max)



#### table for haploid case (with threshold)

number_of_species_table_Haploid_threshold1 <- read.table("number_of_species haploid threshold seed_0.csv", header=FALSE)  
number_of_species_table_Haploid_threshold2 <- read.table("number_of_species haploid threshold seed_2.csv", header=FALSE)
number_of_species_table_Haploid_threshold3 <- read.table("number_of_species haploid threshold seed_3.csv", header=FALSE)  
number_of_species_table_Haploid_threshold4 <- read.table("number_of_species haploid threshold seed_4.csv", header=FALSE)
number_of_species_table_Haploid_threshold5 <- read.table("number_of_species haploid threshold seed_5.csv", header=FALSE)  
number_of_species_table_Haploid_threshold6 <- read.table("number_of_species haploid threshold seed_6.csv", header=FALSE)
number_of_species_table_Haploid_threshold7 <- read.table("number_of_species haploid threshold seed_7.csv", header=FALSE)  
number_of_species_table_Haploid_threshold8 <- read.table("number_of_species haploid threshold seed_8.csv", header=FALSE)
number_of_species_table_Haploid_threshold9 <- read.table("number_of_species haploid threshold seed_9.csv", header=FALSE)  
number_of_species_table_Haploid_threshold10 <- read.table("number_of_species haploid threshold seed_10.csv", header=FALSE)

number_of_species_table_Haploid_threshold <- cbind(number_of_species_table_Haploid_threshold1,
                                                   number_of_species_table_Haploid_threshold2,
                                                   number_of_species_table_Haploid_threshold3,
                                                   number_of_species_table_Haploid_threshold4,
                                                   number_of_species_table_Haploid_threshold5,
                                                   number_of_species_table_Haploid_threshold6,
                                                   number_of_species_table_Haploid_threshold7,
                                                   number_of_species_table_Haploid_threshold8,
                                                   number_of_species_table_Haploid_threshold9,
                                                   number_of_species_table_Haploid_threshold10)

number_of_species_table_Haploid_threshold[1:5,]

(num_of_generations_checked__Haploid_threshold <- dim(number_of_species_table_Haploid_threshold)[1])

# vector of checked generations
checked_generations_vector_Haploid_threshold <- seq(999, 4999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_Haploid_threshold <- apply(X=number_of_species_table_Haploid_threshold, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_Haploid_threshold <- apply(X=number_of_species_table_Haploid_threshold, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_Haploid_threshold <- apply(X=number_of_species_table_Haploid_threshold, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_Haploid_threshold <- apply(X=number_of_species_table_Haploid_threshold, MARGIN=1, FUN=max)



#### table for haploid case with threshold and assortative mating

number_of_species_table_Haploid_threshold_assortative_mating1 <- read.table("number_of_species hap thresh assortative_mating seed_0.csv", header=FALSE)  
number_of_species_table_Haploid_threshold_assortative_mating2 <- read.table("number_of_species hap thresh assortative_mating seed_2.csv", header=FALSE)
number_of_species_table_Haploid_threshold_assortative_mating3 <- read.table("number_of_species hap thresh assortative_mating seed_3.csv", header=FALSE)  
number_of_species_table_Haploid_threshold_assortative_mating4 <- read.table("number_of_species hap thresh assortative_mating seed_4.csv", header=FALSE)
number_of_species_table_Haploid_threshold_assortative_mating5 <- read.table("number_of_species hap thresh assortative_mating seed_5.csv", header=FALSE)  
number_of_species_table_Haploid_threshold_assortative_mating6 <- read.table("number_of_species hap thresh assortative_mating seed_6.csv", header=FALSE)
number_of_species_table_Haploid_threshold_assortative_mating7 <- read.table("number_of_species hap thresh assortative_mating seed_7.csv", header=FALSE)  
number_of_species_table_Haploid_threshold_assortative_mating8 <- read.table("number_of_species hap thresh assortative_mating seed_8.csv", header=FALSE)
number_of_species_table_Haploid_threshold_assortative_mating9 <- read.table("number_of_species hap thresh assortative_mating seed_9.csv", header=FALSE)  
number_of_species_table_Haploid_threshold_assortative_mating10 <- read.table("number_of_species hap thresh assortative_mating seed_10.csv", header=FALSE)

number_of_species_table_Haploid_threshold_assortative_mating <- cbind(number_of_species_table_Haploid_threshold_assortative_mating1,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating2,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating3,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating4,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating5,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating6,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating7,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating8,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating9,
                                                                      number_of_species_table_Haploid_threshold_assortative_mating10)

number_of_species_table_Haploid_threshold_assortative_mating[1:5,]

(num_of_generations_checked__Haploid_threshold_assortative_mating <- dim(number_of_species_table_Haploid_threshold_assortative_mating)[1])

# vector of checked generations
checked_generations_vector_Haploid_threshold_assortative_mating <- seq(999, 4999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_Haploid_threshold_assortative_mating <- apply(X=number_of_species_table_Haploid_threshold_assortative_mating, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_Haploid_threshold_assortative_mating <- apply(X=number_of_species_table_Haploid_threshold_assortative_mating, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_Haploid_threshold_assortative_mating <- apply(X=number_of_species_table_Haploid_threshold_assortative_mating, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_Haploid_threshold_assortative_mating <- apply(X=number_of_species_table_Haploid_threshold_assortative_mating, MARGIN=1, FUN=max)



#### table for sympatric case with assortative mating

number_of_species_table_sympatric_assortative_mating1 <- read.table("number_of_species sympatric assortative mating seed_0.csv", header=FALSE)  
number_of_species_table_sympatric_assortative_mating2 <- read.table("number_of_species sympatric assortative mating seed_2.csv", header=FALSE)
number_of_species_table_sympatric_assortative_mating3 <- read.table("number_of_species sympatric assortative mating seed_3.csv", header=FALSE)  
number_of_species_table_sympatric_assortative_mating4 <- read.table("number_of_species sympatric assortative mating seed_4.csv", header=FALSE)
number_of_species_table_sympatric_assortative_mating5 <- read.table("number_of_species sympatric assortative mating seed_5.csv", header=FALSE)  
number_of_species_table_sympatric_assortative_mating6 <- read.table("number_of_species sympatric assortative mating seed_6.csv", header=FALSE)
number_of_species_table_sympatric_assortative_mating7 <- read.table("number_of_species sympatric assortative mating seed_7.csv", header=FALSE)  
number_of_species_table_sympatric_assortative_mating8 <- read.table("number_of_species sympatric assortative mating seed_8.csv", header=FALSE)
number_of_species_table_sympatric_assortative_mating9 <- read.table("number_of_species sympatric assortative mating seed_9.csv", header=FALSE)  
number_of_species_table_sympatric_assortative_mating10 <- read.table("number_of_species sympatric assortative mating seed_10.csv", header=FALSE)

number_of_species_table_sympatric_assortative_mating <- cbind(number_of_species_table_sympatric_assortative_mating1,
                                                                      number_of_species_table_sympatric_assortative_mating2,
                                                                      number_of_species_table_sympatric_assortative_mating3,
                                                                      number_of_species_table_sympatric_assortative_mating4,
                                                                      number_of_species_table_sympatric_assortative_mating5,
                                                                      number_of_species_table_sympatric_assortative_mating6,
                                                                      number_of_species_table_sympatric_assortative_mating7,
                                                                      number_of_species_table_sympatric_assortative_mating8,
                                                                      number_of_species_table_sympatric_assortative_mating9,
                                                                      number_of_species_table_sympatric_assortative_mating10)

number_of_species_table_sympatric_assortative_mating[1:5,]

(num_of_generations_checked__sympatric_assortative_mating <- dim(number_of_species_table_sympatric_assortative_mating)[1])

# vector of checked generations
checked_generations_vector_sympatric_assortative_mating <- seq(999, 98999, 1000)


# vector of the mean number of species for each checked generation
mean_num_of_species_vector_sympatric_assortative_mating <- apply(X=number_of_species_table_sympatric_assortative_mating, MARGIN=1, FUN=mean)

# vector of the median number of species for each checked generation
median_num_of_species_vector_sympatric_assortative_mating <- apply(X=number_of_species_table_sympatric_assortative_mating, MARGIN=1, FUN=median)

# vector of the min number of species for each checked generation
lower_num_of_species_vector_sympatric_assortative_mating <- apply(X=number_of_species_table_sympatric_assortative_mating, MARGIN=1, FUN=min)

# vector of the max number of species for each checked generation
upper_num_of_species_vector_sympatric_assortative_mating <- apply(X=number_of_species_table_sympatric_assortative_mating, MARGIN=1, FUN=max)






######## Read the number of genetic differences between a random couple of individuals,
######## each from a different species, in generations 5000, 99000 and 149000 in each
######## of the cases for which these generations have been checked.
########
######## Each vector is for for a specific generation in a specific case, and contains
######## the relevant runs (the runs where there was speciation).



#### BDM case


## generation 99000

number_of_differences_table_BDM1_99000 <- read.csv("number_of_differences BDM 99000 seed_0.csv", header=TRUE)  
number_of_differences_table_BDM2_99000 <- read.csv("number_of_differences BDM 99000 seed_2.csv", header=TRUE)  
number_of_differences_table_BDM3_99000 <- read.csv("number_of_differences BDM 99000 seed_3.csv", header=TRUE)  
number_of_differences_table_BDM4_99000 <- read.csv("number_of_differences BDM 99000 seed_4.csv", header=TRUE)  
number_of_differences_table_BDM5_99000 <- read.csv("number_of_differences BDM 99000 seed_5.csv", header=TRUE)  
number_of_differences_table_BDM6_99000 <- read.csv("number_of_differences BDM 99000 seed_6.csv", header=TRUE)  
number_of_differences_table_BDM7_99000 <- read.csv("number_of_differences BDM 99000 seed_7.csv", header=TRUE)  
number_of_differences_table_BDM8_99000 <- read.csv("number_of_differences BDM 99000 seed_8.csv", header=TRUE)  
number_of_differences_table_BDM9_99000 <- read.csv("number_of_differences BDM 99000 seed_9.csv", header=TRUE)  
number_of_differences_table_BDM10_99000 <- read.csv("number_of_differences BDM 99000 seed_10.csv", header=TRUE)  

number_of_differences_vector_BDM_99000 <- c(number_of_differences_table_BDM1_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM2_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM3_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM4_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM5_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM6_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM7_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM8_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM9_99000[1, "number.of.differences"],
                                           number_of_differences_table_BDM10_99000[1, "number.of.differences"])

(number_of_differences_vector_BDM_99000_no_zeroes <- number_of_differences_vector_BDM_99000[number_of_differences_vector_BDM_99000!=0])


## generation 5000

number_of_differences_table_BDM1_5000 <- read.csv("number_of_differences BDM 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_BDM2_5000 <- read.csv("number_of_differences BDM 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_BDM3_5000 <- read.csv("number_of_differences BDM 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_BDM4_5000 <- read.csv("number_of_differences BDM 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_BDM5_5000 <- read.csv("number_of_differences BDM 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_BDM6_5000 <- read.csv("number_of_differences BDM 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_BDM7_5000 <- read.csv("number_of_differences BDM 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_BDM8_5000 <- read.csv("number_of_differences BDM 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_BDM9_5000 <- read.csv("number_of_differences BDM 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_BDM10_5000 <- read.csv("number_of_differences BDM 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_BDM_5000 <- c(number_of_differences_table_BDM1_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM2_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM3_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM4_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM5_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM6_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM7_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM8_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM9_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM10_5000[1, "number.of.differences"])

(number_of_differences_vector_BDM_5000_no_zeroes <- number_of_differences_vector_BDM_5000[number_of_differences_vector_BDM_5000!=0])



#### BDM assortative_mating case


## generation 99000

number_of_differences_table_BDM_assortative_mating1_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_0.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating2_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_2.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating3_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_3.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating4_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_4.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating5_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_5.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating6_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_6.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating7_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_7.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating8_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_8.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating9_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_9.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating10_99000 <- read.csv("number_of_differences BDM_assortative_mating 99000 seed_10.csv", header=TRUE)  

number_of_differences_vector_BDM_assortative_mating_99000 <- c(number_of_differences_table_BDM_assortative_mating1_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating2_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating3_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating4_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating5_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating6_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating7_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating8_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating9_99000[1, "number.of.differences"],
                                            number_of_differences_table_BDM_assortative_mating10_99000[1, "number.of.differences"])

(number_of_differences_vector_BDM_assortative_mating_99000_no_zeroes <- number_of_differences_vector_BDM_assortative_mating_99000[number_of_differences_vector_BDM_assortative_mating_99000!=0])


## generation 5000

number_of_differences_table_BDM_assortative_mating1_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating2_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating3_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating4_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating5_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating6_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating7_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating8_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating9_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_BDM_assortative_mating10_5000 <- read.csv("number_of_differences BDM_assortative_mating 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_BDM_assortative_mating_5000 <- c(number_of_differences_table_BDM_assortative_mating1_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating2_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating3_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating4_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating5_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating6_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating7_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating8_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating9_5000[1, "number.of.differences"],
                                           number_of_differences_table_BDM_assortative_mating10_5000[1, "number.of.differences"])

(number_of_differences_vector_BDM_assortative_mating_5000_no_zeroes <- number_of_differences_vector_BDM_assortative_mating_5000[number_of_differences_vector_BDM_assortative_mating_5000!=0])



#### Diploid threshold case


## generation 149000

number_of_differences_table_Diploid_threshold1_149000 <- read.csv("number_of_differences threshold 149000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold2_149000 <- read.csv("number_of_differences threshold 149000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold3_149000 <- read.csv("number_of_differences threshold 149000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold4_149000 <- read.csv("number_of_differences threshold 149000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold5_149000 <- read.csv("number_of_differences threshold 149000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold6_149000 <- read.csv("number_of_differences threshold 149000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold7_149000 <- read.csv("number_of_differences threshold 149000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold8_149000 <- read.csv("number_of_differences threshold 149000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold9_149000 <- read.csv("number_of_differences threshold 149000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold10_149000 <- read.csv("number_of_differences threshold 149000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_149000 <- c(number_of_differences_table_Diploid_threshold1_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold2_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold3_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold4_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold5_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold6_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold7_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold8_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold9_149000[1, "number.of.differences"],
                                                          number_of_differences_table_Diploid_threshold10_149000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_149000_no_zeroes <- number_of_differences_vector_Diploid_threshold_149000[number_of_differences_vector_Diploid_threshold_149000!=0])


## generation 99000

number_of_differences_table_Diploid_threshold1_99000 <- read.csv("number_of_differences threshold 99000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold2_99000 <- read.csv("number_of_differences threshold 99000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold3_99000 <- read.csv("number_of_differences threshold 99000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold4_99000 <- read.csv("number_of_differences threshold 99000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold5_99000 <- read.csv("number_of_differences threshold 99000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold6_99000 <- read.csv("number_of_differences threshold 99000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold7_99000 <- read.csv("number_of_differences threshold 99000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold8_99000 <- read.csv("number_of_differences threshold 99000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold9_99000 <- read.csv("number_of_differences threshold 99000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold10_99000 <- read.csv("number_of_differences threshold 99000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_99000 <- c(number_of_differences_table_Diploid_threshold1_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold2_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold3_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold4_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold5_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold6_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold7_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold8_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold9_99000[1, "number.of.differences"],
                                            number_of_differences_table_Diploid_threshold10_99000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_99000_no_zeroes <- number_of_differences_vector_Diploid_threshold_99000[number_of_differences_vector_Diploid_threshold_99000!=0])


## generation 5000

number_of_differences_table_Diploid_threshold1_5000 <- read.csv("number_of_differences threshold 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold2_5000 <- read.csv("number_of_differences threshold 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold3_5000 <- read.csv("number_of_differences threshold 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold4_5000 <- read.csv("number_of_differences threshold 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold5_5000 <- read.csv("number_of_differences threshold 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold6_5000 <- read.csv("number_of_differences threshold 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold7_5000 <- read.csv("number_of_differences threshold 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold8_5000 <- read.csv("number_of_differences threshold 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold9_5000 <- read.csv("number_of_differences threshold 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold10_5000 <- read.csv("number_of_differences threshold 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_5000 <- c(number_of_differences_table_Diploid_threshold1_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold2_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold3_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold4_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold5_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold6_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold7_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold8_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold9_5000[1, "number.of.differences"],
                                           number_of_differences_table_Diploid_threshold10_5000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_5000_no_zeroes <- number_of_differences_vector_Diploid_threshold_5000[number_of_differences_vector_Diploid_threshold_5000!=0])



#### Diploid threshold assortative_mating case


## generation 149000

number_of_differences_table_Diploid_threshold_assortative_mating1_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating2_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating3_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating4_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating5_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating6_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating7_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating8_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating9_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating10_149000 <- read.csv("number_of_differences threshold assortative_mating 149000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_assortative_mating_149000 <- c(number_of_differences_table_Diploid_threshold_assortative_mating1_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating2_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating3_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating4_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating5_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating6_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating7_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating8_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating9_149000[1, "number.of.differences"],
                                                                             number_of_differences_table_Diploid_threshold_assortative_mating10_149000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_assortative_mating_149000_no_zeroes <- number_of_differences_vector_Diploid_threshold_assortative_mating_149000[number_of_differences_vector_Diploid_threshold_assortative_mating_149000!=0])


## generation 99000

number_of_differences_table_Diploid_threshold_assortative_mating1_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating2_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating3_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating4_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating5_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating6_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating7_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating8_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating9_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating10_99000 <- read.csv("number_of_differences threshold assortative_mating 99000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_assortative_mating_99000 <- c(number_of_differences_table_Diploid_threshold_assortative_mating1_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating2_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating3_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating4_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating5_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating6_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating7_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating8_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating9_99000[1, "number.of.differences"],
                                                               number_of_differences_table_Diploid_threshold_assortative_mating10_99000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_assortative_mating_99000_no_zeroes <- number_of_differences_vector_Diploid_threshold_assortative_mating_99000[number_of_differences_vector_Diploid_threshold_assortative_mating_99000!=0])


## generation 5000

number_of_differences_table_Diploid_threshold_assortative_mating1_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating2_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating3_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating4_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating5_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating6_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating7_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating8_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating9_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_Diploid_threshold_assortative_mating10_5000 <- read.csv("number_of_differences threshold assortative_mating 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Diploid_threshold_assortative_mating_5000 <- c(number_of_differences_table_Diploid_threshold_assortative_mating1_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating2_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating3_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating4_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating5_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating6_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating7_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating8_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating9_5000[1, "number.of.differences"],
                                                              number_of_differences_table_Diploid_threshold_assortative_mating10_5000[1, "number.of.differences"])

(number_of_differences_vector_Diploid_threshold_assortative_mating_5000_no_zeroes <- number_of_differences_vector_Diploid_threshold_assortative_mating_5000[number_of_differences_vector_Diploid_threshold_assortative_mating_5000!=0])



#### Haploid threshold case


## generation 5000

number_of_differences_table_Haploid_threshold1_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold2_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold3_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold4_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold5_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold6_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold7_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold8_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold9_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold10_5000 <- read.csv("number_of_differences haploid threshold 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Haploid_threshold_5000 <- c(number_of_differences_table_Haploid_threshold1_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold2_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold3_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold4_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold5_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold6_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold7_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold8_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold9_5000[1, "number.of.differences"],
                                                         number_of_differences_table_Haploid_threshold10_5000[1, "number.of.differences"])

(number_of_differences_vector_Haploid_threshold_5000_no_zeroes <- number_of_differences_vector_Haploid_threshold_5000[number_of_differences_vector_Haploid_threshold_5000!=0])



#### Haploid threshold assortative_mating case


## generation 5000

number_of_differences_table_Haploid_threshold_assortative_mating1_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_0.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating2_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_2.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating3_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_3.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating4_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_4.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating5_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_5.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating6_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_6.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating7_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_7.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating8_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_8.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating9_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_9.csv", header=TRUE)  
number_of_differences_table_Haploid_threshold_assortative_mating10_5000 <- read.csv("number_of_differences haploid threshold assortative_mating 5000 seed_10.csv", header=TRUE)  

number_of_differences_vector_Haploid_threshold_assortative_mating_5000 <- c(number_of_differences_table_Haploid_threshold_assortative_mating1_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating2_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating3_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating4_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating5_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating6_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating7_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating8_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating9_5000[1, "number.of.differences"],
                                                                            number_of_differences_table_Haploid_threshold_assortative_mating10_5000[1, "number.of.differences"])

(number_of_differences_vector_Haploid_threshold_assortative_mating_5000_no_zeroes <- number_of_differences_vector_Haploid_threshold_assortative_mating_5000[number_of_differences_vector_Haploid_threshold_assortative_mating_5000!=0])





######## Read tables of division into species in each checked generation for the cases
######## of BDM and BDM_assortative_mating.



#### table for the regular BDM case

division_into_species_table_BDM1 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_0.csv", header=TRUE)  
division_into_species_table_BDM2 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_2.csv", header=TRUE)
division_into_species_table_BDM3 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_3.csv", header=TRUE)  
division_into_species_table_BDM4 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_4.csv", header=TRUE)
division_into_species_table_BDM5 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_5.csv", header=TRUE)  
division_into_species_table_BDM6 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_6.csv", header=TRUE)
division_into_species_table_BDM7 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_7.csv", header=TRUE)  
division_into_species_table_BDM8 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_8.csv", header=TRUE)
division_into_species_table_BDM9 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_9.csv", header=TRUE)  
division_into_species_table_BDM10 <- read.csv("circle BDM table_for_species 2014_December_10 mut_01 theta_02 ngen_150000 seed_10.csv", header=TRUE)


#remove headers

colnames(division_into_species_table_BDM1) <- NULL
colnames(division_into_species_table_BDM2) <- NULL
colnames(division_into_species_table_BDM3) <- NULL
colnames(division_into_species_table_BDM4) <- NULL
colnames(division_into_species_table_BDM5) <- NULL
colnames(division_into_species_table_BDM6) <- NULL
colnames(division_into_species_table_BDM7) <- NULL
colnames(division_into_species_table_BDM8) <- NULL
colnames(division_into_species_table_BDM9) <- NULL
colnames(division_into_species_table_BDM10) <- NULL

division_into_species_table_BDM1 <- as.matrix(division_into_species_table_BDM1)
division_into_species_table_BDM2 <- as.matrix(division_into_species_table_BDM2)
division_into_species_table_BDM3 <- as.matrix(division_into_species_table_BDM3)  
division_into_species_table_BDM4 <- as.matrix(division_into_species_table_BDM4)
division_into_species_table_BDM5 <- as.matrix(division_into_species_table_BDM5)  
division_into_species_table_BDM6 <- as.matrix(division_into_species_table_BDM6)
division_into_species_table_BDM7 <- as.matrix(division_into_species_table_BDM7)  
division_into_species_table_BDM8 <- as.matrix(division_into_species_table_BDM8)
division_into_species_table_BDM9 <- as.matrix(division_into_species_table_BDM9)
division_into_species_table_BDM10 <- as.matrix(division_into_species_table_BDM10)

max_num_of_rows <- max(dim(division_into_species_table_BDM1)[1], dim(division_into_species_table_BDM2)[1], 
                       dim(division_into_species_table_BDM3)[1], dim(division_into_species_table_BDM4)[1], 
                       dim(division_into_species_table_BDM5)[1], dim(division_into_species_table_BDM6)[1],
                       dim(division_into_species_table_BDM7)[1], dim(division_into_species_table_BDM8)[1], 
                       dim(division_into_species_table_BDM9)[1], dim(division_into_species_table_BDM10)[1])

if (max_num_of_rows > dim(division_into_species_table_BDM1)[1])
   division_into_species_table_BDM1 <- rbind(division_into_species_table_BDM1, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM1)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM1)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM2)[1])
   division_into_species_table_BDM2 <- rbind(division_into_species_table_BDM2, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM2)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM2)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM3)[1])
   division_into_species_table_BDM3 <- rbind(division_into_species_table_BDM3, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM3)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM3)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM4)[1])
   division_into_species_table_BDM4 <- rbind(division_into_species_table_BDM4, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM4)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM4)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM5)[1])
   division_into_species_table_BDM5 <- rbind(division_into_species_table_BDM5, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM5)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM5)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM6)[1])
   division_into_species_table_BDM6 <- rbind(division_into_species_table_BDM6, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM6)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM6)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM7)[1])
   division_into_species_table_BDM7 <- rbind(division_into_species_table_BDM7, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM7)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM7)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM8)[1])
   division_into_species_table_BDM8 <- rbind(division_into_species_table_BDM8, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM8)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM8)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM9)[1])
   division_into_species_table_BDM9 <- rbind(division_into_species_table_BDM9, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM9)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM9)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM10)[1])
   division_into_species_table_BDM10 <- rbind(division_into_species_table_BDM10, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM10)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM10)[1])))


# table of division into species in each checked generation and for each seed (run)
# for the case of BDM.
# Each run appears on the right of the former.
division_into_species_table_BDM <- cbind(division_into_species_table_BDM1,
                                         division_into_species_table_BDM2,
                                         division_into_species_table_BDM3,
                                         division_into_species_table_BDM4,
                                         division_into_species_table_BDM5,
                                         division_into_species_table_BDM6,
                                         division_into_species_table_BDM7,
                                         division_into_species_table_BDM8,
                                         division_into_species_table_BDM9,
                                         division_into_species_table_BDM10)

dim(division_into_species_table_BDM)
#  2503   30 

nrows_original <- dim(division_into_species_table_BDM)[1]

division_into_species_table_BDM[1, ]

# the max number of species in all generations and runs
(max_num_of_species_BDM <- max(number_of_species_table_BDM))
# 7



#### Clean the table from NA's and put the runs one below the other.
####
#### Add a column of run number and a column of generation number (in thousands), and remove
#### the first column which is redundant.
####
#### Add a column of species number in {1, 2, 3, ..., max_num_of_species_BDM} to
#### each of the species in each generation.
#### In each generation, different numbers are attached to different species.
####
#### Start by cleaning each run's table separately.

# Matrix of first row of each generation.
# Each row is for a different generation.
# Each column is for a different run.
first_row_of_each_generation <- matrix(rep(0, num_of_generations_checked__BDM * sample_size), ncol=sample_size)


## Run 1

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM1[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM1 <- cbind(rep(1, nrows_original), division_into_species_table_BDM1, rep(0, nrows_original))
colnames(division_into_species_table_BDM1) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 1
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


division_into_species_table_BDM1 <- data.frame(run = as.numeric(division_into_species_table_BDM1[, "run"]),
                    gen = as.numeric(division_into_species_table_BDM1[, "gen"]),
                    section.size = as.numeric(division_into_species_table_BDM1[, "section size"]),
                    name.of.species = as.character(division_into_species_table_BDM1[, "name of species"]),
                    species.number = as.numeric(division_into_species_table_BDM1[, "species number"]))
                    


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM1[current_row, "species.number"]))){
  
  
  # if the last row was an NA row
  if (NA_row){
    
    # if the current row is an NA row
    if (is.na(division_into_species_table_BDM1[current_row, "section.size"])){
      
      # delete the current row
      division_into_species_table_BDM1 <- division_into_species_table_BDM1[-current_row, ]
      
    } else{
      # if the current row contains information  
      
      current_gen <- current_gen + 1
      first_row_of_each_generation[current_gen, current_run] <- current_row
      division_into_species_table_BDM1[current_row, "species.number"] <- 1
      division_into_species_table_BDM1[current_row, "gen"] <- current_gen
      current_species_number <- 1
      current_row <- current_row + 1
      NA_row <- FALSE
      
    }
    
  } else{
    # if the last row contained information
    
    # if the current row is an NA row
    if (is.na(division_into_species_table_BDM1[current_row, "section.size"])){
      number_of_species_in_current_generation <- current_species_number
      if (current_gen > 1){
        
        # Find a permutation of the species numbers in the current generation
        # that minimizes the number of individuals who changed species relative
        # to the previous generation.
        permutation_for_minimal_change <- 1:number_of_species_in_current_generation
        minimal_change <- 0
        current_individual_to_check <- 1
        next_individual_to_check <- 1
        current_section_previous_generation <- 1
        current_section_current_generation <- 1
        
        backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
        backwards_found <- FALSE
        while (!backwards_found){
          if(is.na(division_into_species_table_BDM1[first_row_of_each_generation[(current_gen - backwards), current_run], "species.number"])){
            backwards <- backwards + 1
          } else{
            backwards_found <- TRUE
          }
        }
        
        current_num_of_species_previous_generation <- 
          division_into_species_table_BDM1[first_row_of_each_generation[(current_gen - backwards), current_run], "species.number"]
        current_num_of_species_current_generation <- 1
        last_individual_in_current_section_previous_generation <-
          division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section.size"]
        last_individual_in_current_section_current_generation <-
          division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run]), "section.size"]
        
        while (current_individual_to_check <= num_of_individuals_checked){
          next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                              as.numeric(last_individual_in_current_section_current_generation))
          if (division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                current_section_previous_generation - 1), "species.number"] !=
              division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen), current_run] + 
                                                current_section_current_generation - 1), "species.number"]){
            minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
          }
          if (last_individual_in_current_section_previous_generation < 
              as.numeric(last_individual_in_current_section_current_generation)){
            current_num_of_species_previous_generation <- 
              division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                  current_section_previous_generation), "species.number"]
            current_section_previous_generation <- current_section_previous_generation + 1
            last_individual_in_current_section_previous_generation <- 
              last_individual_in_current_section_previous_generation + 
              division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                  current_section_previous_generation - 1), "section.size"]
          } else{
            if (last_individual_in_current_section_previous_generation >
                as.numeric(last_individual_in_current_section_current_generation)){
              current_num_of_species_current_generation <- 
                division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] +
                                                    current_section_current_generation), "species.number"]
              current_section_current_generation <- current_section_current_generation + 1
              last_individual_in_current_section_current_generation <- 
                last_individual_in_current_section_current_generation + 
                division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                    current_section_current_generation - 1), "section.size"]
            } else{
              current_num_of_species_previous_generation <- 
                division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                    current_section_previous_generation), "species.number"]
              current_section_previous_generation <- current_section_previous_generation + 1
              last_individual_in_current_section_previous_generation <- 
                last_individual_in_current_section_previous_generation + 
                division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                    current_section_previous_generation - 1), "section.size"]
              current_num_of_species_current_generation <- 
                division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                    current_section_current_generation), "species.number"]
              current_section_current_generation <- current_section_current_generation + 1
              last_individual_in_current_section_current_generation <- 
                last_individual_in_current_section_current_generation + 
                division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                    current_section_current_generation - 1), "section.size"]
            }
          }
          current_individual_to_check <- next_individual_to_check
        }
        if (number_of_species_in_current_generation > 1){
          if (number_of_species_in_current_generation < 9){
            possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
            if (number_of_species_in_current_generation > 2){
              nrow_possible_permutations <- dim(possible_permutations)[1]
            } else{
              nrow_possible_permutations <- 1
              possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
            }
            for (permutation_number in 1:nrow_possible_permutations){
              change_in_this_permutation <- 0
              current_individual_to_check <- 1
              next_individual_to_check <- 1
              current_section_previous_generation <- 1
              current_section_current_generation <- 1
              last_individual_in_current_section_previous_generation <- division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section.size"]
              last_individual_in_current_section_current_generation <- division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run]), "section.size"]
              current_num_of_species_previous_generation <- division_into_species_table_BDM1[first_row_of_each_generation[(current_gen - backwards), current_run], "species.number"]
              current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
              while (current_individual_to_check <= num_of_individuals_checked){
                next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                as.numeric(last_individual_in_current_section_current_generation)) + 1
                if (division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species.number"] !=
                    possible_permutations[permutation_number, division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species.number"]]){
                  change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                }
                if (last_individual_in_current_section_previous_generation < 
                    as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                    division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                        current_section_previous_generation), "species.number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                    last_individual_in_current_section_previous_generation + 
                    division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                        current_section_previous_generation - 1), "section.size"]
                } else{
                  if (last_individual_in_current_section_previous_generation >
                      as.numeric(last_individual_in_current_section_current_generation)){
                    current_num_of_species_current_generation <- 
                      division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] +
                                                          current_section_current_generation), "species.number"]
                    current_section_current_generation <- current_section_current_generation + 1
                    last_individual_in_current_section_current_generation <- 
                      last_individual_in_current_section_current_generation + 
                      division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                          current_section_current_generation - 1), "section.size"]
                  } else{
                    current_num_of_species_previous_generation <- 
                      division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                          current_section_previous_generation), "species.number"]
                    current_section_previous_generation <- current_section_previous_generation + 1
                    last_individual_in_current_section_previous_generation <- 
                      last_individual_in_current_section_previous_generation + 
                      division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                          current_section_previous_generation - 1), "section.size"]
                    current_num_of_species_current_generation <- 
                      division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                          current_section_current_generation), "species.number"]
                    current_section_current_generation <- current_section_current_generation + 1
                    last_individual_in_current_section_current_generation <- 
                      last_individual_in_current_section_current_generation + 
                      division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                          current_section_current_generation - 1), "section.size"]
                  }
                }
                current_individual_to_check <- next_individual_to_check
              } 
              if (change_in_this_permutation < minimal_change){
                minimal_change <- change_in_this_permutation
                permutation_for_minimal_change <- possible_permutations[permutation_number, ]
              }
            }
            
            for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
              division_into_species_table_BDM1[j, "species.number"] <- permutation_for_minimal_change[division_into_species_table_BDM1[j, "species.number"]]
            } 
          } else{
            # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
            
            for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
              division_into_species_table_BDM1[j, "species.number"] <- NA
            }
          }   
        }
      }
      # delete the current row
      division_into_species_table_BDM1 <- division_into_species_table_BDM1[-current_row, ]
      NA_row <- TRUE
      
    } else{
      # if the current row contains information 
      
      species_match_found <- FALSE
      current_row_to_check <- current_row - 1
      
      while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM1[current_row_to_check, "gen"] == current_gen)){
        if (division_into_species_table_BDM1[current_row, "name.of.species"] ==
            division_into_species_table_BDM1[current_row_to_check, "name.of.species"]){
          division_into_species_table_BDM1[current_row, "species.number"] <-
            division_into_species_table_BDM1[current_row_to_check, "species.number"]
          species_match_found <- TRUE
        } else{
          current_row_to_check <- current_row_to_check - 1
        }
      }           
      
      if (!species_match_found){
        current_species_number <- current_species_number + 1
        division_into_species_table_BDM1[current_row, "species.number"] <- current_species_number
      }
      division_into_species_table_BDM1[current_row, "gen"] <- current_gen
      current_row <- current_row + 1
      NA_row <- FALSE
      
    }  
    
  }
  
}



































while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM1[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM1[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM1 <- division_into_species_table_BDM1[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM1[current_row, "species number"] <- 1
         division_into_species_table_BDM1[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM1[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM1
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               as.numeric(division_into_species_table_BDM1[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
              as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"])
            last_individual_in_current_section_previous_generation <- as.numeric(last_individual_in_current_section_previous_generation)
            last_individual_in_current_section_current_generation <-
              as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run]), "section size"])
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(as.numeric(last_individual_in_current_section_current_generation)))
               if (division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"])
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"])
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"])
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"])
                  } else{
                     current_num_of_species_previous_generation <- 
                        as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"])
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                       as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"])
                     current_num_of_species_current_generation <- 
                       as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"])
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                       as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"])
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"])
                     last_individual_in_current_section_current_generation <- as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run]), "section size"])
                     current_num_of_species_previous_generation <- as.numeric(division_into_species_table_BDM1[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"]) !=
                               possible_permutations[permutation_number, as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"])]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                             as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"])
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                             as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"])
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"])
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"])
                           } else{
                              current_num_of_species_previous_generation <- 
                                as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"])
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"])
                              current_num_of_species_current_generation <- 
                                as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"])
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 as.numeric(division_into_species_table_BDM1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"])
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM1[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM1[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM1[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM1 <- division_into_species_table_BDM1[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM1[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM1[current_row, "name of species"] ==
                   division_into_species_table_BDM1[current_row_to_check, "name of species"]){
               division_into_species_table_BDM1[current_row, "species number"] <-
                  division_into_species_table_BDM1[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM1[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM1[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 2

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM2[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM2 <- cbind(rep(2, nrows_original), division_into_species_table_BDM2, rep(0, nrows_original))
colnames(division_into_species_table_BDM2) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 2
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM2[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM2[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM2 <- division_into_species_table_BDM2[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM2[current_row, "species number"] <- 1
         division_into_species_table_BDM2[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM2[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM2
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM2[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM2[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM2[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM2[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM2[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM2 <- division_into_species_table_BDM2[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM2[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM2[current_row, "name of species"] ==
                   division_into_species_table_BDM2[current_row_to_check, "name of species"]){
               division_into_species_table_BDM2[current_row, "species number"] <-
                  division_into_species_table_BDM2[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM2[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM2[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 3

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM3[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM3 <- cbind(rep(3, nrows_original), division_into_species_table_BDM3, rep(0, nrows_original))
colnames(division_into_species_table_BDM3) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 3
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM3[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM3[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM3 <- division_into_species_table_BDM3[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM3[current_row, "species number"] <- 1
         division_into_species_table_BDM3[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM3[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM3
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM3[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM3[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM3[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM3[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM3[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM3 <- division_into_species_table_BDM3[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM3[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM3[current_row, "name of species"] ==
                   division_into_species_table_BDM3[current_row_to_check, "name of species"]){
               division_into_species_table_BDM3[current_row, "species number"] <-
                  division_into_species_table_BDM3[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM3[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM3[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 4

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM4[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM4 <- cbind(rep(4, nrows_original), division_into_species_table_BDM4, rep(0, nrows_original))
colnames(division_into_species_table_BDM4) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 4
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM4[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM4[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM4 <- division_into_species_table_BDM4[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM4[current_row, "species number"] <- 1
         division_into_species_table_BDM4[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM4[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM4
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM4[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM4[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM4[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM4[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM4[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM4 <- division_into_species_table_BDM4[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM4[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM4[current_row, "name of species"] ==
                   division_into_species_table_BDM4[current_row_to_check, "name of species"]){
               division_into_species_table_BDM4[current_row, "species number"] <-
                  division_into_species_table_BDM4[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM4[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM4[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 5

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM5[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM5 <- cbind(rep(5, nrows_original), division_into_species_table_BDM5, rep(0, nrows_original))
colnames(division_into_species_table_BDM5) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 5
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM5[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM5[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM5 <- division_into_species_table_BDM5[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM5[current_row, "species number"] <- 1
         division_into_species_table_BDM5[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM5[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM5
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM5[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM5[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM5[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM5[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM5[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM5 <- division_into_species_table_BDM5[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM5[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM5[current_row, "name of species"] ==
                   division_into_species_table_BDM5[current_row_to_check, "name of species"]){
               division_into_species_table_BDM5[current_row, "species number"] <-
                  division_into_species_table_BDM5[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM5[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM5[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 6

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM6[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM6 <- cbind(rep(6, nrows_original), division_into_species_table_BDM6, rep(0, nrows_original))
colnames(division_into_species_table_BDM6) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 6
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM6[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM6[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM6 <- division_into_species_table_BDM6[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM6[current_row, "species number"] <- 1
         division_into_species_table_BDM6[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM6[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM6
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM6[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM6[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM6[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM6[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM6[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM6 <- division_into_species_table_BDM6[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM6[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM6[current_row, "name of species"] ==
                   division_into_species_table_BDM6[current_row_to_check, "name of species"]){
               division_into_species_table_BDM6[current_row, "species number"] <-
                  division_into_species_table_BDM6[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM6[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM6[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 7

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM7[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM7 <- cbind(rep(7, nrows_original), division_into_species_table_BDM7, rep(0, nrows_original))
colnames(division_into_species_table_BDM7) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 7
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM7[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM7[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM7 <- division_into_species_table_BDM7[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM7[current_row, "species number"] <- 1
         division_into_species_table_BDM7[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM7[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM7
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM7[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM7[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM7[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM7[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM7[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM7 <- division_into_species_table_BDM7[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM7[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM7[current_row, "name of species"] ==
                   division_into_species_table_BDM7[current_row_to_check, "name of species"]){
               division_into_species_table_BDM7[current_row, "species number"] <-
                  division_into_species_table_BDM7[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM7[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM7[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 8

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM8[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM8 <- cbind(rep(8, nrows_original), division_into_species_table_BDM8, rep(0, nrows_original))
colnames(division_into_species_table_BDM8) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 8
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM8[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM8[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM8 <- division_into_species_table_BDM8[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM8[current_row, "species number"] <- 1
         division_into_species_table_BDM8[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM8[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM8
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM8[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (as.numeric(last_individual_in_current_section_previous_generation) < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     as.character(as.numeric(last_individual_in_current_section_previous_generation) + 
                     as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]))
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        as.character(as.numeric(last_individual_in_current_section_current_generation) + 
                        as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]))
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        as.character(as.numeric(last_individual_in_current_section_previous_generation) + 
                        as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]))
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        as.character(as.numeric(last_individual_in_current_section_current_generation) + 
                        as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]))
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM8[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              as.character(as.numeric(last_individual_in_current_section_previous_generation) + 
                              as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]))
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 as.character(as.numeric(last_individual_in_current_section_current_generation) + 
                                 as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]))
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 as.character(as.numeric(last_individual_in_current_section_previous_generation) + 
                                 as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]))
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 as.character(as.numeric(last_individual_in_current_section_current_generation) + 
                                 as.numeric(division_into_species_table_BDM8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]))
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM8[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM8[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM8[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM8 <- division_into_species_table_BDM8[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM8[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM8[current_row, "name of species"] ==
                   division_into_species_table_BDM8[current_row_to_check, "name of species"]){
               division_into_species_table_BDM8[current_row, "species number"] <-
                  division_into_species_table_BDM8[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM8[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM8[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 9

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM9[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM9 <- cbind(rep(9, nrows_original), division_into_species_table_BDM9, rep(0, nrows_original))
colnames(division_into_species_table_BDM9) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 9
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM9[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM9[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM9 <- division_into_species_table_BDM9[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM9[current_row, "species number"] <- 1
         division_into_species_table_BDM9[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM9[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM9
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM9[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM9[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM9[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM9[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM9[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM9 <- division_into_species_table_BDM9[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM9[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM9[current_row, "name of species"] ==
                   division_into_species_table_BDM9[current_row_to_check, "name of species"]){
               division_into_species_table_BDM9[current_row, "species number"] <-
                  division_into_species_table_BDM9[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM9[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM9[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 10

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM10[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM10 <- cbind(rep(10, nrows_original), division_into_species_table_BDM10, rep(0, nrows_original))
colnames(division_into_species_table_BDM10) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 10
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM) && !((current_gen == num_of_generations_checked__BDM) && is.na(division_into_species_table_BDM10[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM10[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM10 <- division_into_species_table_BDM10[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM10[current_row, "species number"] <- 1
         division_into_species_table_BDM10[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM10[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM10
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM10[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                            current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                               current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                              current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                              current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] +
                                                                                 current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                 current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                 current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                 current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                 current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                 current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM10[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                       current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                       current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] +
                                                                                          current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                          current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                          current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                          current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                          current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                          current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM10[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM10[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM10[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM10 <- division_into_species_table_BDM10[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM10[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM10[current_row, "name of species"] ==
                   division_into_species_table_BDM10[current_row_to_check, "name of species"]){
               division_into_species_table_BDM10[current_row, "species number"] <-
                  division_into_species_table_BDM10[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM10[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM10[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}



#### Full table of division into species in each checked generation and for each seed (run)
#### for the case of BDM.
#### Each run appears below the former.

division_into_species_table_BDM <- rbind(division_into_species_table_BDM1,
                                         division_into_species_table_BDM2,
                                         division_into_species_table_BDM3,
                                         division_into_species_table_BDM4,
                                         division_into_species_table_BDM5,
                                         division_into_species_table_BDM6,
                                         division_into_species_table_BDM7,
                                         division_into_species_table_BDM8,
                                         division_into_species_table_BDM9,
                                         division_into_species_table_BDM10)



#### Write table of division into species to Excel

write.csv(division_into_species_table_BDM, "paper 2015/tables for species/division_into_species_table_BDM.csv", row.names=FALSE)



#### table for the regular BDM_assortative_mating case

division_into_species_table_BDM_assortative_mating1 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_0.csv", header=FALSE)  
division_into_species_table_BDM_assortative_mating2 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_2.csv", header=FALSE)
division_into_species_table_BDM_assortative_mating3 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_3.csv", header=FALSE)  
division_into_species_table_BDM_assortative_mating4 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_4.csv", header=FALSE)
division_into_species_table_BDM_assortative_mating5 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_5.csv", header=FALSE)  
division_into_species_table_BDM_assortative_mating6 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_6.csv", header=FALSE)
division_into_species_table_BDM_assortative_mating7 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_7.csv", header=FALSE)  
division_into_species_table_BDM_assortative_mating8 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_8.csv", header=FALSE)
division_into_species_table_BDM_assortative_mating9 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_9.csv", header=FALSE)  
division_into_species_table_BDM_assortative_mating10 <- read.csv("paper 2015/tables for species/circle assortative_mating table_for_species 2014_December_10 mut_01 theta_02 ngen_100000 seed_10.csv", header=FALSE)



max_num_of_rows <- max(dim(division_into_species_table_BDM_assortative_mating1)[1], dim(division_into_species_table_BDM_assortative_mating2)[1], 
                       dim(division_into_species_table_BDM_assortative_mating3)[1], dim(division_into_species_table_BDM_assortative_mating4)[1], 
                       dim(division_into_species_table_BDM_assortative_mating5)[1], dim(division_into_species_table_BDM_assortative_mating6)[1],
                       dim(division_into_species_table_BDM_assortative_mating7)[1], dim(division_into_species_table_BDM_assortative_mating8)[1], 
                       dim(division_into_species_table_BDM_assortative_mating9)[1], dim(division_into_species_table_BDM_assortative_mating10)[1])

if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating1)[1])
  division_into_species_table_BDM_assortative_mating1 <- rbind(division_into_species_table_BDM_assortative_mating1, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating1)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating1)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating2)[1])
  division_into_species_table_BDM_assortative_mating2 <- rbind(division_into_species_table_BDM_assortative_mating2, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating2)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating2)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating3)[1])
  division_into_species_table_BDM_assortative_mating3 <- rbind(division_into_species_table_BDM_assortative_mating3, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating3)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating3)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating4)[1])
  division_into_species_table_BDM_assortative_mating4 <- rbind(division_into_species_table_BDM_assortative_mating4, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating4)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating4)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating5)[1])
  division_into_species_table_BDM_assortative_mating5 <- rbind(division_into_species_table_BDM_assortative_mating5, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating5)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating5)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating6)[1])
  division_into_species_table_BDM_assortative_mating6 <- rbind(division_into_species_table_BDM_assortative_mating6, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating6)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating6)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating7)[1])
  division_into_species_table_BDM_assortative_mating7 <- rbind(division_into_species_table_BDM_assortative_mating7, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating7)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating7)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating8)[1])
  division_into_species_table_BDM_assortative_mating8 <- rbind(division_into_species_table_BDM_assortative_mating8, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating8)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating8)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating9)[1])
  division_into_species_table_BDM_assortative_mating9 <- rbind(division_into_species_table_BDM_assortative_mating9, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating9)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating9)[1])))
if (max_num_of_rows > dim(division_into_species_table_BDM_assortative_mating10)[1])
  division_into_species_table_BDM_assortative_mating10 <- rbind(division_into_species_table_BDM_assortative_mating10, matrix(rep(NA, 3*(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating10)[1])), nrow=(max_num_of_rows - dim(division_into_species_table_BDM_assortative_mating10)[1])))


# table of division into species in each checked generation and for each seed (run)
# for the case of BDM_assortative_mating.
# Each run appears on the right of the former.
division_into_species_table_BDM_assortative_mating <- cbind(division_into_species_table_BDM_assortative_mating1,
                                         division_into_species_table_BDM_assortative_mating2,
                                         division_into_species_table_BDM_assortative_mating3,
                                         division_into_species_table_BDM_assortative_mating4,
                                         division_into_species_table_BDM_assortative_mating5,
                                         division_into_species_table_BDM_assortative_mating6,
                                         division_into_species_table_BDM_assortative_mating7,
                                         division_into_species_table_BDM_assortative_mating8,
                                         division_into_species_table_BDM_assortative_mating9,
                                         division_into_species_table_BDM_assortative_mating10)

dim(division_into_species_table_BDM_assortative_mating)
#  2128   30 

nrows_original <- dim(division_into_species_table_BDM_assortative_mating)[1]

division_into_species_table_BDM_assortative_mating[1, ]

# the max number of species in all generations and runs
(max_num_of_species_BDM_assortative_mating <- max(number_of_species_table_BDM_assortative_mating))
# 13



#### Clean the table from NA's and put the runs one below the other.
####
#### Add a column of run number and a column of generation number (in thousands), and remove
#### the first column which is redundant.
####
#### Add a column of species number in {1, 2, 3, ..., max_num_of_species_BDM_assortative_mating} to
#### each of the species in each generation.
#### In each generation, different numbers are attached to different species.
####
#### Start by cleaning each run's table separately.

# Matrix of first row of each generation.
# Each row is for a different generation.
# Each column is for a different run.
first_row_of_each_generation <- matrix(rep(0, num_of_generations_checked__BDM_assortative_mating * sample_size), ncol=sample_size)


## Run 1

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating1[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating1 <- cbind(rep(1, nrows_original), division_into_species_table_BDM_assortative_mating1, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating1) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 1
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating1[current_row, "species number"]))){
  
   # if the last row was an NA row
   if (NA_row){
    
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating1[current_row, "section size"])){
      
         # delete the current row
         division_into_species_table_BDM_assortative_mating1 <- division_into_species_table_BDM_assortative_mating1[-current_row, ]
      
      } else{
      # if the current row contains information  
      
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating1[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating1[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
      
      }
    
   } else{
   # if the last row contained information
    
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating1[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
        
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
                    
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating1
                          [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
        
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating1[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                        
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                      as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                     current_section_previous_generation - 1), "species number"] !=
                        division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen), current_run] + 
                           current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                     as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                        current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                        current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                        as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] +
                           current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                           current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                           current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                           division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                              current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                           current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                           division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                              current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating1[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                            possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating1[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
              
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating1[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating1[j, "species number"]]
                  } 
               } else{
               # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating1[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating1 <- division_into_species_table_BDM_assortative_mating1[-current_row, ]
         NA_row <- TRUE
      
      } else{
      # if the current row contains information 
      
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
      
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating1[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating1[current_row, "name of species"] ==
                  division_into_species_table_BDM_assortative_mating1[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating1[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating1[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
      
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating1[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating1[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
      
      }  
    
   }
  
}


## Run 2

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating2[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating2 <- cbind(rep(2, nrows_original), division_into_species_table_BDM_assortative_mating2, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating2) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 2
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating2[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating2[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating2 <- division_into_species_table_BDM_assortative_mating2[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating2[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating2[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating2[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating2
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating2[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating2[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating2[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating2[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating2[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating2[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating2 <- division_into_species_table_BDM_assortative_mating2[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating2[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating2[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating2[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating2[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating2[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating2[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating2[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 3

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating3[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating3 <- cbind(rep(3, nrows_original), division_into_species_table_BDM_assortative_mating3, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating3) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 3
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating3[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating3[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating3 <- division_into_species_table_BDM_assortative_mating3[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating3[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating3[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating3[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating3
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating3[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating3[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating3[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating3[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating3[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating3[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating3 <- division_into_species_table_BDM_assortative_mating3[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating3[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating3[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating3[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating3[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating3[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating3[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating3[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 4

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating4[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating4 <- cbind(rep(4, nrows_original), division_into_species_table_BDM_assortative_mating4, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating4) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 4
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating4[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating4[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating4 <- division_into_species_table_BDM_assortative_mating4[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating4[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating4[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating4[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating4
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating4[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating4[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating4[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating4[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating4[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating4[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating4 <- division_into_species_table_BDM_assortative_mating4[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating4[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating4[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating4[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating4[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating4[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating4[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating4[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 5

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating5[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating5 <- cbind(rep(5, nrows_original), division_into_species_table_BDM_assortative_mating5, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating5) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 5
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating5[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating5[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating5 <- division_into_species_table_BDM_assortative_mating5[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating5[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating5[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating5[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating5
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating5[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating5[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating5[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating5[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating5[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating5[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating5 <- division_into_species_table_BDM_assortative_mating5[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating5[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating5[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating5[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating5[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating5[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating5[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating5[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 6

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating6[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating6 <- cbind(rep(6, nrows_original), division_into_species_table_BDM_assortative_mating6, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating6) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 6
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating6[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating6[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating6 <- division_into_species_table_BDM_assortative_mating6[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating6[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating6[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating6[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating6
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating6[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating6[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating6[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating6[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating6[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating6[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating6 <- division_into_species_table_BDM_assortative_mating6[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating6[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating6[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating6[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating6[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating6[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating6[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating6[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 7

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating7[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating7 <- cbind(rep(7, nrows_original), division_into_species_table_BDM_assortative_mating7, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating7) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 7
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating7[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating7[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating7 <- division_into_species_table_BDM_assortative_mating7[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating7[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating7[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating7[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating7
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating7[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating7[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating7[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating7[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating7[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating7[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating7 <- division_into_species_table_BDM_assortative_mating7[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating7[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating7[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating7[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating7[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating7[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating7[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating7[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 8

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating8[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating8 <- cbind(rep(8, nrows_original), division_into_species_table_BDM_assortative_mating8, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating8) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 8
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating8[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating8[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating8 <- division_into_species_table_BDM_assortative_mating8[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating8[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating8[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating8[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating8
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating8[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating8[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating8[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating8[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating8[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating8[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating8 <- division_into_species_table_BDM_assortative_mating8[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating8[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating8[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating8[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating8[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating8[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating8[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating8[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 9

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating9[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating9 <- cbind(rep(9, nrows_original), division_into_species_table_BDM_assortative_mating9, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating9) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 9
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating9[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating9[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating9 <- division_into_species_table_BDM_assortative_mating9[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating9[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating9[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating9[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating9
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating9[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating9[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating9[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating9[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating9[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating9[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating9 <- division_into_species_table_BDM_assortative_mating9[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating9[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating9[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating9[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating9[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating9[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating9[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating9[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}


## Run 10

# Instead of the first (redundant) column, put a column of run number and a
# column of generation number (in thousands).
# Add a column of species number.
division_into_species_table_BDM_assortative_mating10[, 1] <- rep(0, nrows_original)
division_into_species_table_BDM_assortative_mating10 <- cbind(rep(10, nrows_original), division_into_species_table_BDM_assortative_mating10, rep(0, nrows_original))
colnames(division_into_species_table_BDM_assortative_mating10) <- c("run", "gen", "section size", "name of species", "species number")

current_run <- 10
current_gen <- 0
current_row <- 1

# a boolean to mark whether the current row is an NA row
NA_row <- TRUE


while ((current_gen <= num_of_generations_checked__BDM_assortative_mating) && !((current_gen == num_of_generations_checked__BDM_assortative_mating) && is.na(division_into_species_table_BDM_assortative_mating10[current_row, "species number"]))){
   
   # if the last row was an NA row
   if (NA_row){
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating10[current_row, "section size"])){
         
         # delete the current row
         division_into_species_table_BDM_assortative_mating10 <- division_into_species_table_BDM_assortative_mating10[-current_row, ]
         
      } else{
         # if the current row contains information  
         
         current_gen <- current_gen + 1
         first_row_of_each_generation[current_gen, current_run] <- current_row
         division_into_species_table_BDM_assortative_mating10[current_row, "species number"] <- 1
         division_into_species_table_BDM_assortative_mating10[current_row, "gen"] <- current_gen
         current_species_number <- 1
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }
      
   } else{
      # if the last row contained information
      
      # if the current row is an NA row
      if (is.na(division_into_species_table_BDM_assortative_mating10[current_row, "section size"])){
         number_of_species_in_current_generation <- current_species_number
         if (current_gen > 1){
            
            # Find a permutation of the species numbers in the current generation
            # that minimizes the number of individuals who changed species relative
            # to the previous generation.
            permutation_for_minimal_change <- 1:number_of_species_in_current_generation
            minimal_change <- 0
            current_individual_to_check <- 1
            next_individual_to_check <- 1
            current_section_previous_generation <- 1
            current_section_current_generation <- 1
            
            backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
            backwards_found <- FALSE
            while (!backwards_found){
               if(is.na(division_into_species_table_BDM_assortative_mating10
                        [first_row_of_each_generation[(current_gen - backwards), current_run], "species number"])){
                  backwards <- backwards + 1
               } else{
                  backwards_found <- TRUE
               }
            }
            
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM_assortative_mating10[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
            current_num_of_species_current_generation <- 1
            last_individual_in_current_section_previous_generation <-
               division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
            last_individual_in_current_section_current_generation <-
               division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run]), "section size"]
            
            while (current_individual_to_check <= num_of_individuals_checked){
               next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                   as.numeric(last_individual_in_current_section_current_generation))
               if (division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                           current_section_previous_generation - 1), "species number"] !=
                      division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen), current_run] + 
                                                                              current_section_current_generation - 1), "species number"]){
                  minimal_change <- minimal_change + (next_individual_to_check - current_individual_to_check)
               }
               if (last_individual_in_current_section_previous_generation < 
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                             current_section_previous_generation - 1), "section size"]
               } else{
                  if (last_individual_in_current_section_previous_generation >
                         as.numeric(last_individual_in_current_section_current_generation)){
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] +
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  } else{
                     current_num_of_species_previous_generation <- 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation), "species number"]
                     current_section_previous_generation <- current_section_previous_generation + 1
                     last_individual_in_current_section_previous_generation <- 
                        last_individual_in_current_section_previous_generation + 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                current_section_previous_generation - 1), "section size"]
                     current_num_of_species_current_generation <- 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation), "species number"]
                     current_section_current_generation <- current_section_current_generation + 1
                     last_individual_in_current_section_current_generation <- 
                        last_individual_in_current_section_current_generation + 
                        division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                current_section_current_generation - 1), "section size"]
                  }
               }
               current_individual_to_check <- next_individual_to_check
            }
            if (number_of_species_in_current_generation > 1){
               if (number_of_species_in_current_generation < 9){
                  possible_permutations <- allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))
                  if (number_of_species_in_current_generation > 2){
                     nrow_possible_permutations <- dim(possible_permutations)[1]
                  } else{
                     nrow_possible_permutations <- 1
                     possible_permutations <- t(as.matrix(allPerms(number_of_species_in_current_generation, control = how(maxperm=99999))))
                  }
                  for (permutation_number in 1:nrow_possible_permutations){
                     change_in_this_permutation <- 0
                     current_individual_to_check <- 1
                     next_individual_to_check <- 1
                     current_section_previous_generation <- 1
                     current_section_current_generation <- 1
                     last_individual_in_current_section_previous_generation <- division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run]), "section size"]
                     last_individual_in_current_section_current_generation <- division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run]), "section size"]
                     current_num_of_species_previous_generation <- division_into_species_table_BDM_assortative_mating10[first_row_of_each_generation[(current_gen - backwards), current_run], "species number"]
                     current_num_of_species_current_generation <- possible_permutations[permutation_number, 1]
                     while (current_individual_to_check <= num_of_individuals_checked){
                        next_individual_to_check <- min(as.numeric(last_individual_in_current_section_previous_generation),
                                                        as.numeric(last_individual_in_current_section_current_generation)) + 1
                        if (division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + current_section_previous_generation - 1), "species number"] !=
                               possible_permutations[permutation_number, division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen), current_run] + current_section_current_generation - 1), "species number"]]){
                           change_in_this_permutation <- change_in_this_permutation + (next_individual_to_check - current_individual_to_check)
                        }
                        if (last_individual_in_current_section_previous_generation < 
                               as.numeric(last_individual_in_current_section_current_generation)){
                           current_num_of_species_previous_generation <- 
                              division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation), "species number"]
                           current_section_previous_generation <- current_section_previous_generation + 1
                           last_individual_in_current_section_previous_generation <- 
                              last_individual_in_current_section_previous_generation + 
                              division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                      current_section_previous_generation - 1), "section size"]
                        } else{
                           if (last_individual_in_current_section_previous_generation >
                                  as.numeric(last_individual_in_current_section_current_generation)){
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] +
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           } else{
                              current_num_of_species_previous_generation <- 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation), "species number"]
                              current_section_previous_generation <- current_section_previous_generation + 1
                              last_individual_in_current_section_previous_generation <- 
                                 last_individual_in_current_section_previous_generation + 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[(current_gen - backwards), current_run] + 
                                                                                         current_section_previous_generation - 1), "section size"]
                              current_num_of_species_current_generation <- 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation), "species number"]
                              current_section_current_generation <- current_section_current_generation + 1
                              last_individual_in_current_section_current_generation <- 
                                 last_individual_in_current_section_current_generation + 
                                 division_into_species_table_BDM_assortative_mating10[(first_row_of_each_generation[current_gen, current_run] + 
                                                                                         current_section_current_generation - 1), "section size"]
                           }
                        }
                        current_individual_to_check <- next_individual_to_check
                     } 
                     if (change_in_this_permutation < minimal_change){
                        minimal_change <- change_in_this_permutation
                        permutation_for_minimal_change <- possible_permutations[permutation_number, ]
                     }
                  }
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating10[j, "species number"] <- permutation_for_minimal_change[division_into_species_table_BDM_assortative_mating10[j, "species number"]]
                  } 
               } else{
                  # The number of species in this generation is >=9, and therefore we cannot check all permutations.   
                  
                  for (j in first_row_of_each_generation[current_gen, current_run]:(current_row - 1)){
                     division_into_species_table_BDM_assortative_mating10[j, "species number"] <- NA
                  }
               }   
            }
         }
         # delete the current row
         division_into_species_table_BDM_assortative_mating10 <- division_into_species_table_BDM_assortative_mating10[-current_row, ]
         NA_row <- TRUE
         
      } else{
         # if the current row contains information 
         
         species_match_found <- FALSE
         current_row_to_check <- current_row - 1
         
         while ((current_row_to_check > 0) && (!species_match_found) && (division_into_species_table_BDM_assortative_mating10[current_row_to_check, "gen"] == current_gen)){
            if (division_into_species_table_BDM_assortative_mating10[current_row, "name of species"] ==
                   division_into_species_table_BDM_assortative_mating10[current_row_to_check, "name of species"]){
               division_into_species_table_BDM_assortative_mating10[current_row, "species number"] <-
                  division_into_species_table_BDM_assortative_mating10[current_row_to_check, "species number"]
               species_match_found <- TRUE
            } else{
               current_row_to_check <- current_row_to_check - 1
            }
         }           
         
         if (!species_match_found){
            current_species_number <- current_species_number + 1
            division_into_species_table_BDM_assortative_mating10[current_row, "species number"] <- current_species_number
         }
         division_into_species_table_BDM_assortative_mating10[current_row, "gen"] <- current_gen
         current_row <- current_row + 1
         NA_row <- FALSE
         
      }  
      
   }
   
}



#### Full table of division into species in each checked generation and for each seed (run)
#### for the case of BDM_assortative_mating.
#### Each run appears below the former.

division_into_species_table_BDM_assortative_mating <- rbind(division_into_species_table_BDM_assortative_mating1,
                                         division_into_species_table_BDM_assortative_mating2,
                                         division_into_species_table_BDM_assortative_mating3,
                                         division_into_species_table_BDM_assortative_mating4,
                                         division_into_species_table_BDM_assortative_mating5,
                                         division_into_species_table_BDM_assortative_mating6,
                                         division_into_species_table_BDM_assortative_mating7,
                                         division_into_species_table_BDM_assortative_mating8,
                                         division_into_species_table_BDM_assortative_mating9,
                                         division_into_species_table_BDM_assortative_mating10)



#### Write table of division into species to Excel

write.csv(division_into_species_table_BDM_assortative_mating, "paper 2015/tables for species/division_into_species_table_BDM_assortative_mating.csv", row.names=FALSE)




######## Perform resampling tests for the absolute value (to make it a two-sided test) of
######## the statistic integral_of_difference_between_medians,  
######## absolute_integral_of_difference_between_medians, to check the null hypothesis
######## that the integral is 0, against the alternative that the number of species in
######## one of the speciation types tends to have larger values than the other. 



#### BDM vs. BDM_assortative_mating


# number of permutations
R = 10000

# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM, number_of_species_table_BDM_assortative_mating)


for (i in 1:R) {
  
   group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
   group1_matrix <- unified_matrix[, group1_indices]
   group2_matrix <- unified_matrix[, -group1_indices]
  
   median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
   median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
   integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
   absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM - median_num_of_species_vector_BDM_assortative_mating))
# -144

(p_value__BDM__BDM_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0.0015
#  =>  Number of species in  BDM_assortative_mating > BDM



#### BDM vs. Diploid_threshold


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM, number_of_species_table_Diploid_threshold[1:59,])


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM - median_num_of_species_vector_Diploid_threshold[1:59]))
# 90

(p_value__BDM__Diploid_threshold <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  BDM > Diploid_threshold



#### BDM vs. Diploid_threshold_assortative_mating vector


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM, number_of_species_table_Diploid_threshold_assortative_mating[1:59,])


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM - median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59]))
# 90

(p_value__BDM__Diploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0.0825
# Cannot conclude whether  Diploid_threshold_assortative_mating != BDM



#### BDM_assortative_mating vs. Diploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM_assortative_mating, number_of_species_table_Diploid_threshold_assortative_mating[1:59,])


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM_assortative_mating - median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59]))
# 234

(p_value__BDM_assortative_mating__Diploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0.0166
#  =>  Number of species in  BDM_assortative_mating > Diploid_threshold_assortative_mating



#### Diploid_threshold vs. Diploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_Diploid_threshold, number_of_species_table_Diploid_threshold_assortative_mating)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_Diploid_threshold - median_num_of_species_vector_Diploid_threshold_assortative_mating))
# -1152

(p_value__Diploid_threshold__Diploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Diploid_threshold_assortative_mating > Diploid_threshold



#### BDM vs. Haploid_threshold vector


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM[1:5,], number_of_species_table_Haploid_threshold)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM[1:5] - median_num_of_species_vector_Haploid_threshold))
# -82.5

(p_value__BDM__Haploid_threshold <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold > BDM



#### Diploid_threshold vs. Haploid_threshold vector


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_Diploid_threshold[1:5,], number_of_species_table_Haploid_threshold)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_Diploid_threshold[1:5] - median_num_of_species_vector_Haploid_threshold))
# -82.5

(p_value__Diploid_threshold__Haploid_threshold <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold > Diploid_threshold



#### Diploid_threshold_assortative_mating vs. Haploid_threshold


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_Diploid_threshold_assortative_mating[1:5,], number_of_species_table_Haploid_threshold)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:5] - median_num_of_species_vector_Haploid_threshold))
# -82.5

(p_value__Diploid_threshold_assortative_mating__Haploid_threshold <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold > Diploid_threshold_assortative_mating



#### BDM vs. Haploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM[1:5,], number_of_species_table_Haploid_threshold_assortative_mating)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM[1:5] - median_num_of_species_vector_Haploid_threshold_assortative_mating))
# -495

(p_value__BDM__Haploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold_assortative_mating > BDM



#### BDM_assortative_mating vs. Haploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_BDM_assortative_mating[1:5,], number_of_species_table_Haploid_threshold_assortative_mating)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_BDM_assortative_mating[1:5] - median_num_of_species_vector_Haploid_threshold_assortative_mating))
# -466.5

(p_value__BDM_assortative_mating__Haploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold_assortative_mating > BDM_assortative_mating



#### Diploid_threshold_assortative_mating vs. Haploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_Diploid_threshold_assortative_mating[1:5,], number_of_species_table_Haploid_threshold_assortative_mating)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:5] - median_num_of_species_vector_Haploid_threshold_assortative_mating))
# -495

(p_value__Diploid_threshold_assortative_mating__Haploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold_assortative_mating > Diploid_threshold_assortative_mating



#### Haploid_threshold vs. Haploid_threshold_assortative_mating


# vector of integral_of_difference_between_medians results
integral_of_difference_between_medians <- rep(0, R)

# vector of absolute_integral_of_difference_between_medians results
absolute_integral_of_difference_between_medians <- rep(0, R)

# unified matrix of number of species for the two groups
unified_matrix <- cbind(number_of_species_table_Haploid_threshold, number_of_species_table_Haploid_threshold_assortative_mating)


for (i in 1:R) {
  
  group1_indices <- sample(1:(2*sample_size), size=sample_size, replace=F)
  
  group1_matrix <- unified_matrix[, group1_indices]
  group2_matrix <- unified_matrix[, -group1_indices]
  
  median_num_of_species_vector_group1 <- apply(X=group1_matrix, MARGIN=1, FUN=median)
  median_num_of_species_vector_group2 <- apply(X=group2_matrix, MARGIN=1, FUN=median)
  
  integral_of_difference_between_medians[i] <- sum(median_num_of_species_vector_group1 - median_num_of_species_vector_group2)
  absolute_integral_of_difference_between_medians[i] <- abs(integral_of_difference_between_medians[i])
  
}

(original_groups_integral_of_difference_between_medians <- sum(median_num_of_species_vector_Haploid_threshold - median_num_of_species_vector_Haploid_threshold_assortative_mating))
# -412.5

(p_value__Haploid_threshold__Haploid_threshold_assortative_mating <- (sum(absolute_integral_of_difference_between_medians  >  abs(original_groups_integral_of_difference_between_medians)) / R))
# p-value = 0
#  =>  Number of species in  Haploid_threshold_assortative_mating > Haploid_threshold










######## plots of number of species over time: medians and range


checked_generations_vector_BDM_0 <- c(0, checked_generations_vector_BDM)

(lower_num_of_species_vector_BDM <- c(1, lower_num_of_species_vector_BDM))
(median_num_of_species_vector_BDM <- c(1, median_num_of_species_vector_BDM))
(upper_num_of_species_vector_BDM <- c(1, upper_num_of_species_vector_BDM))

(lower_num_of_species_vector_BDM_assortative_mating <- c(1, lower_num_of_species_vector_BDM_assortative_mating))
(median_num_of_species_vector_BDM_assortative_mating <- c(1, median_num_of_species_vector_BDM_assortative_mating))
(upper_num_of_species_vector_BDM_assortative_mating <- c(1, upper_num_of_species_vector_BDM_assortative_mating))

(lower_num_of_species_vector_Diploid_threshold <- c(1, lower_num_of_species_vector_Diploid_threshold))
(median_num_of_species_vector_Diploid_threshold <- c(1, median_num_of_species_vector_Diploid_threshold))
(upper_num_of_species_vector_Diploid_threshold <- c(1, upper_num_of_species_vector_Diploid_threshold))

(lower_num_of_species_vector_Diploid_threshold_assortative_mating <- c(1, lower_num_of_species_vector_Diploid_threshold_assortative_mating))
(median_num_of_species_vector_Diploid_threshold_assortative_mating <- c(1, median_num_of_species_vector_Diploid_threshold_assortative_mating))
(upper_num_of_species_vector_Diploid_threshold_assortative_mating <- c(1, upper_num_of_species_vector_Diploid_threshold_assortative_mating))

(lower_num_of_species_vector_Haploid_threshold <- c(1, lower_num_of_species_vector_Haploid_threshold))
(median_num_of_species_vector_Haploid_threshold <- c(1, median_num_of_species_vector_Haploid_threshold))
(upper_num_of_species_vector_Haploid_threshold <- c(1, upper_num_of_species_vector_Haploid_threshold))

(lower_num_of_species_vector_Haploid_threshold_assortative_mating <- c(1, lower_num_of_species_vector_Haploid_threshold_assortative_mating))
(median_num_of_species_vector_Haploid_threshold_assortative_mating <- c(1, median_num_of_species_vector_Haploid_threshold_assortative_mating))
(upper_num_of_species_vector_Haploid_threshold_assortative_mating <- c(1, upper_num_of_species_vector_Haploid_threshold_assortative_mating))

(lower_num_of_species_vector_sympatric_assortative_mating <- c(1, lower_num_of_species_vector_sympatric_assortative_mating))
(median_num_of_species_vector_sympatric_assortative_mating <- c(1, median_num_of_species_vector_sympatric_assortative_mating))
(upper_num_of_species_vector_sympatric_assortative_mating <- c(1, upper_num_of_species_vector_sympatric_assortative_mating))




## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
# (plot Diploid log range.png, plot Diploid log range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
lo <- loess((log(c(median_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(c(median_num_of_species_vector_BDM) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(upper_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(lower_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_BDM_50 <- c(0,checked_generations_vector_BDM + 50)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM_50, y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(upper_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(lower_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_100 <- c(0,checked_generations_vector_Diploid_threshold + 100)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(median_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(upper_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(lower_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"), col=c("black", "blue", "red", "purple"), fill=c("black", "blue", "red", "purple"))





## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
y_lim_min <- 0.5
y_lim_max <- 1 + max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating)))
exp(y_lim_max - 1)

y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)

# Set custom x-axis ticks
x_ticks <- c(0, 20000, 40000, 60000, 80000, 100000)
x_labels <- format(x_ticks, big.mark = ",", scientific = FALSE)

# Helper: lower bound for log(1 + 1)
log_min <- log(1 + 1)

# BDM
lo <- loess((log(c(median_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span = 0.4)
plot(x = checked_generations_vector_BDM_0, y = pmax(predict(lo), log_min),
     xlab = 'generation', ylab = 'number of species',
     main = 'BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating',
     cex.main = 0.8, xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)
axis(2, at = y_tic_marks_location, labels = y_tic_marks_values, las = 2)
axis(1, at = x_ticks, labels = x_labels)

lo <- loess((log(c(upper_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_BDM_0, y = pmax(predict(lo), log_min), ann = FALSE,
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

lo <- loess((log(c(lower_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_BDM_0, y = pmax(predict(lo), log_min), ann = FALSE,
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# BDM assortative mating
checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ checked_generations_vector_BDM_50, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_BDM_50, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lwd = 2)

lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ checked_generations_vector_BDM_50, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_BDM_50, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ checked_generations_vector_BDM_50, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_BDM_50, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# Diploid threshold
checked_generations_vector_Diploid_threshold_100 <- c(0, checked_generations_vector_Diploid_threshold + 100)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold + 1)) ~ checked_generations_vector_Diploid_threshold_100, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_100, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lwd = 2)

lo <- loess((log(upper_num_of_species_vector_Diploid_threshold + 1)) ~ checked_generations_vector_Diploid_threshold_100, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_100, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

lo <- loess((log(lower_num_of_species_vector_Diploid_threshold + 1)) ~ checked_generations_vector_Diploid_threshold_100, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_100, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# Diploid threshold assortative mating
checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ checked_generations_vector_Diploid_threshold_150, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_150, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lwd = 2)

lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ checked_generations_vector_Diploid_threshold_150, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_150, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ checked_generations_vector_Diploid_threshold_150, span = 0.4)
par(new = TRUE)
plot(x = checked_generations_vector_Diploid_threshold_150, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 90000), ylim = c(y_lim_min, y_lim_max), yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

legend("topleft",
       legend = c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"),
       col = c("black", "blue", "red", "purple"),
       fill = c("black", "blue", "red", "purple"))





## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
y_lim_min <- 0.5
y_lim_max <- 1 + max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating)))
exp(y_lim_max - 1)

y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)

# Set custom x-axis ticks
x_ticks <- c(0, 20000, 40000, 60000, 80000, 100000)
x_labels <- format(x_ticks, big.mark = ",", scientific = FALSE)

# Helper: lower bound for log(1 + 1)
log_min <- log(1 + 1)

# Truncation helper
truncate <- function(x, y, max_x = 99009) {
  idx <- which(x <= max_x)
  list(x = x[idx], y = y[idx])
}

# BDM
tmp <- truncate(checked_generations_vector_BDM_0, median_num_of_species_vector_BDM)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
plot(x = tmp$x, y = pmax(predict(lo), log_min),
     xlab = 'generation', ylab = 'number of species',
     main = 'BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating',
     cex.main = 0.8, xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)
axis(2, at = y_tic_marks_location, labels = y_tic_marks_values, las = 2)
axis(1, at = x_ticks, labels = x_labels)

tmp <- truncate(checked_generations_vector_BDM_0, upper_num_of_species_vector_BDM)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE,
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

tmp <- truncate(checked_generations_vector_BDM_0, lower_num_of_species_vector_BDM)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE,
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# BDM assortative mating
checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)

tmp <- truncate(checked_generations_vector_BDM_50, median_num_of_species_vector_BDM_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)

tmp <- truncate(checked_generations_vector_BDM_50, upper_num_of_species_vector_BDM_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

tmp <- truncate(checked_generations_vector_BDM_50, lower_num_of_species_vector_BDM_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'blue',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# Diploid threshold
checked_generations_vector_Diploid_threshold_100 <- c(0, checked_generations_vector_Diploid_threshold + 100)

tmp <- truncate(checked_generations_vector_Diploid_threshold_100, median_num_of_species_vector_Diploid_threshold)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)

tmp <- truncate(checked_generations_vector_Diploid_threshold_100, upper_num_of_species_vector_Diploid_threshold)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

tmp <- truncate(checked_generations_vector_Diploid_threshold_100, lower_num_of_species_vector_Diploid_threshold)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'red',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

# Diploid threshold assortative mating
checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)

tmp <- truncate(checked_generations_vector_Diploid_threshold_150, median_num_of_species_vector_Diploid_threshold_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)

tmp <- truncate(checked_generations_vector_Diploid_threshold_150, upper_num_of_species_vector_Diploid_threshold_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

tmp <- truncate(checked_generations_vector_Diploid_threshold_150, lower_num_of_species_vector_Diploid_threshold_assortative_mating)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
par(new = TRUE)
plot(x = tmp$x, y = pmax(predict(lo), log_min), ann = FALSE, col = 'purple',
     xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lty = "longdash")

legend("topleft",
       legend = c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"),
       col = c("black", "blue", "red", "purple"),
       fill = c("black", "blue", "red", "purple"))




## The Final Version
## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
y_lim_min <- 0.5
y_lim_max <- 1 + max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating)))
exp(y_lim_max - 1)

y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)

# Set custom x-axis ticks
x_ticks <- c(0, 20000, 40000, 60000, 80000, 100000)
x_labels <- format(x_ticks, big.mark = ",", scientific = FALSE)

# Lower bound for y-axis in log space
log_min <- log(1 + 1)

# Truncation helper
truncate <- function(x, y, max_x = 99009) {
  idx <- which(x <= max_x)
  list(x = x[idx], y = y[idx])
}

# Plot loess curve that starts smoothly at y = log(2)
plot_loess_with_smooth_start <- function(x, y, span = 0.4, col = "black", lty = 1, lwd = 2, new = TRUE) {
  tmp <- truncate(x, y)
  lo <- loess(log(tmp$y + 1) ~ tmp$x, span = span)
  predicted <- predict(lo)
  
  # Smooth vertical shift so curve starts at log(2)
  offset <- log_min - predicted[1]
  adjusted <- pmax(predicted + offset, log_min)
  
  if (new) par(new = TRUE)
  lines(x = tmp$x, y = adjusted, col = col, lty = lty, lwd = lwd)
}

# BDM Median – First plot
tmp <- truncate(checked_generations_vector_BDM_0, median_num_of_species_vector_BDM)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
predicted <- predict(lo)
offset <- log_min - predicted[1]
adjusted <- pmax(predicted + offset, log_min)

plot(x = tmp$x, y = adjusted,
     xlab = 'generation', ylab = 'number of species',
     main = 'BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating',
     cex.main = 0.8, xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 2)
axis(2, at = y_tic_marks_location, labels = y_tic_marks_values, las = 2)
axis(1, at = x_ticks, labels = x_labels)

# BDM upper/lower
plot_loess_with_smooth_start(checked_generations_vector_BDM_0, upper_num_of_species_vector_BDM, lty = "longdash")
plot_loess_with_smooth_start(checked_generations_vector_BDM_0, lower_num_of_species_vector_BDM, lty = "longdash")

# BDM assortative mating
checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, median_num_of_species_vector_BDM_assortative_mating, col = "blue", lwd = 2)
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, upper_num_of_species_vector_BDM_assortative_mating, col = "blue", lty = "longdash")
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, lower_num_of_species_vector_BDM_assortative_mating, col = "blue", lty = "longdash")

# Diploid threshold
checked_generations_vector_Diploid_threshold_100 <- c(0, checked_generations_vector_Diploid_threshold + 100)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, median_num_of_species_vector_Diploid_threshold, col = "red", lwd = 2)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, upper_num_of_species_vector_Diploid_threshold, col = "red", lty = "longdash")
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, lower_num_of_species_vector_Diploid_threshold, col = "red", lty = "longdash")

# Diploid threshold assortative mating
checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, median_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lwd = 2)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, upper_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lty = "longdash")
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, lower_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lty = "longdash")

# Legend
legend("topleft",
       legend = c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"),
       col = c("black", "blue", "red", "purple"),
       fill = c("black", "blue", "red", "purple"))





## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
y_lim_min <- 0.5
y_lim_max <- 1 + max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating)))
exp(y_lim_max - 1)

y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)

# Set custom x-axis ticks
x_ticks <- c(0, 20000, 40000, 60000, 80000, 100000)
x_labels <- format(x_ticks, big.mark = ",", scientific = FALSE)

# Lower bound for y-axis in log space
log_min <- log(1 + 1)

# Truncation helper
truncate <- function(x, y, max_x = 99009) {
  idx <- which(x <= max_x)
  list(x = x[idx], y = y[idx])
}

# Plot loess curve that starts smoothly at y = log(2)
plot_loess_with_smooth_start <- function(x, y, span = 0.4, col = "black", lty = 1, lwd = 1, new = TRUE) {
  tmp <- truncate(x, y)
  lo <- loess(log(tmp$y + 1) ~ tmp$x, span = span)
  predicted <- predict(lo)
  
  # Smooth vertical shift so curve starts at log(2)
  offset <- log_min - predicted[1]
  adjusted <- pmax(predicted + offset, log_min)
  
  if (new) par(new = TRUE)
  lines(x = tmp$x, y = adjusted, col = col, lty = lty, lwd = lwd)
}

# BDM Median – First plot
tmp <- truncate(checked_generations_vector_BDM_0, median_num_of_species_vector_BDM)
lo <- loess(log(tmp$y + 1) ~ tmp$x, span = 0.4)
predicted <- predict(lo)
offset <- log_min - predicted[1]
adjusted <- pmax(predicted + offset, log_min)

plot(x = tmp$x, y = adjusted,
     xlab = 'generation', ylab = 'number of species',
     main = 'BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating',
     cex.main = 0.8, xlim = c(0, 100000), ylim = c(y_lim_min, y_lim_max),
     yaxt = "n", xaxt = "n", type = "l", lwd = 1)
axis(2, at = y_tic_marks_location, labels = y_tic_marks_values, las = 2)
axis(1, at = x_ticks, labels = x_labels)

# BDM upper/lower
plot_loess_with_smooth_start(checked_generations_vector_BDM_0, upper_num_of_species_vector_BDM, lty = "longdash", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_BDM_0, lower_num_of_species_vector_BDM, lty = "longdash", lwd = 1)

# BDM assortative mating
checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, median_num_of_species_vector_BDM_assortative_mating, col = "blue", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, upper_num_of_species_vector_BDM_assortative_mating, col = "blue", lty = "longdash", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_BDM_50, lower_num_of_species_vector_BDM_assortative_mating, col = "blue", lty = "longdash", lwd = 1)

# Diploid threshold
checked_generations_vector_Diploid_threshold_100 <- c(0, checked_generations_vector_Diploid_threshold + 100)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, median_num_of_species_vector_Diploid_threshold, col = "red", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, upper_num_of_species_vector_Diploid_threshold, col = "red", lty = "longdash", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_100, lower_num_of_species_vector_Diploid_threshold, col = "red", lty = "longdash", lwd = 1)

# Diploid threshold assortative mating
checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, median_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, upper_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lty = "longdash", lwd = 1)
plot_loess_with_smooth_start(checked_generations_vector_Diploid_threshold_150, lower_num_of_species_vector_Diploid_threshold_assortative_mating, col = "purple", lty = "longdash", lwd = 1)

# Legend
legend("topleft",
       legend = c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"),
       col = c("black", "blue", "red", "purple"),
       fill = c("black", "blue", "red", "purple"))










## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating

# log scale on y axis with range
# (plot Diploid log range.png, plot Diploid log range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
lo <- loess((log(c(median_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(c(median_num_of_species_vector_BDM) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(upper_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(lower_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_BDM_50 <- c(0,checked_generations_vector_BDM + 50)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM_50, y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(upper_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(lower_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_100 <- c(0,checked_generations_vector_Diploid_threshold + 100)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(median_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(upper_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(lower_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_150 <- c(0, checked_generations_vector_Diploid_threshold + 150)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 99000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating"), col=c("black", "blue", "red", "purple"), fill=c("black", "blue", "red", "purple"))









## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating 99000 generations

# log scale on y axis with range
# (plot Diploid log range 99000.png, plot Diploid log range 99000.eps)
# plot_Diploid_log_range_99000.TIFF, plot_Diploid_log_range_99000.eps, plot_Diploid_log_range_99000.pdf
median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- median_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
lower_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
median_num_of_species_vector_Diploid_threshold_99000 <- median_num_of_species_vector_Diploid_threshold[1:100]
upper_num_of_species_vector_Diploid_threshold_99000 <- upper_num_of_species_vector_Diploid_threshold[1:100]
lower_num_of_species_vector_Diploid_threshold_99000 <- lower_num_of_species_vector_Diploid_threshold[1:100]

checked_generations_vector_Diploid_threshold_99000 <- c(0,checked_generations_vector_Diploid_threshold[1:99])

y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
lo <- loess((log(c(median_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=pmax(predict(lo),1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(c(median_num_of_species_vector_BDM) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(upper_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(lower_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM_50, y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(upper_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(lower_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_99000_100 <- checked_generations_vector_Diploid_threshold_99000 + 100
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(median_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=log(upper_num_of_species_vector_Diploid_threshold_99000 + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=log(lower_num_of_species_vector_Diploid_threshold_99000 + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_99000_150 <- checked_generations_vector_Diploid_threshold_99000 + 150
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
par(new=T)
##plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
par(new=T)
##plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("BDM (1)", "BDM assortative mating (2)", "Diploid threshold (3)", "Diploid threshold assortative mating (4)"), col=c("black", "blue", "red", "purple"), fill=c("black", "blue", "red", "purple"))



## BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating 99000 generations

# log scale on y axis with range
# (plot Diploid log range 99000.png, plot Diploid log range 99000.eps)
# plot_Diploid_log_range_99000.TIFF, plot_Diploid_log_range_99000.eps, plot_Diploid_log_range_99000.pdf
median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- median_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
lower_num_of_species_vector_Diploid_threshold_assortative_mating_99000 <- lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:100]
median_num_of_species_vector_Diploid_threshold_99000 <- median_num_of_species_vector_Diploid_threshold[1:100]
upper_num_of_species_vector_Diploid_threshold_99000 <- upper_num_of_species_vector_Diploid_threshold[1:100]
lower_num_of_species_vector_Diploid_threshold_99000 <- lower_num_of_species_vector_Diploid_threshold[1:100]

checked_generations_vector_Diploid_threshold_99000 <- c(0,checked_generations_vector_Diploid_threshold[1:99])

y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
lo <- loess((log(c(median_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=pmax(predict(lo),1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(c(median_num_of_species_vector_BDM) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(upper_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM) + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold, y=log(c(lower_num_of_species_vector_BDM, rep(NA, (149-59))) + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold & Diploid threshold and assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM_50, y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(upper_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_BDM_50, y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_50, y=log(c(lower_num_of_species_vector_BDM_assortative_mating, rep(NA, (149-59))) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_99000_100 <- checked_generations_vector_Diploid_threshold_99000 + 100
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(median_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
par(new=T)
##plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=log(upper_num_of_species_vector_Diploid_threshold_99000 + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_99000 + 1)) ~ (checked_generations_vector_Diploid_threshold_99000_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_99000_100, y=log(lower_num_of_species_vector_Diploid_threshold_99000 + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_99000_150 <- checked_generations_vector_Diploid_threshold_99000 + 150
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
par(new=T)
##plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating_99000 + 1)) ~ ((checked_generations_vector_Diploid_threshold_99000_150)), span=0.4)
par(new=T)
##plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_99000_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("BDM (1)", "BDM assortative mating (2)", "Diploid threshold (3)", "Diploid threshold assortative mating (4)"), col=c("black", "blue", "red", "purple"), fill=c("black", "blue", "red", "purple"))





## Diploid threshold &  Diploid threshold assortative_mating

# log scale on y axis with range
# (plot Diploid with threshold log range.png, plot Diploid with threshold log range.eps)
# plot_diploid_with_threshold_log_range.pdf, plot_diploid_with_threshold_log_range.TIFF, plot_diploid_with_threshold_log_range.eps
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Diploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
checked_generations_vector_Diploid_threshold_100 <- c(0,checked_generations_vector_Diploid_threshold + 100)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), xlab='generation', ylab='number of species', main='Diploid threshold & Diploid threshold and assortative mating', cex.main=0.85, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(median_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(upper_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold + 1)) ~ (checked_generations_vector_Diploid_threshold_100), span=0.4)
par(new=T)
plot(x=checked_generations_vector_Diploid_threshold_100, y=predict(lo), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Diploid_threshold_100, y=log(lower_num_of_species_vector_Diploid_threshold + 1), ann=F, col='red', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Diploid_threshold_150 <- c(0,checked_generations_vector_Diploid_threshold + 150)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Diploid_threshold_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_Diploid_threshold_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Diploid_threshold_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating + 1), ann=F, col='purple', xlim=c(0, 150000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("Diploid threshold (3)", "Diploid threshold assortative mating (4)"), col=c("red", "purple"), fill=c("red", "purple"))


## BDM, Diploid threshold, Haploid threshold

# log scale on y axis with range
# (plot Haploid threshold and BDM log range.png, plot Haploid threshold and BDM log range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Haploid_threshold))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50)
lo <- loess((log(median_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='BDM, Diploid threshold & Haploid threshold', cex.main=0.85, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(median_num_of_species_vector_BDM + 1), xlab='generation', ylab='number of species', main='BDM, Diploid threshold & Haploid threshold', cex.main=0.85, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(upper_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(lower_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_BDM_50 <- c(0,checked_generations_vector_BDM + 50)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_50)), span=0.4)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_50), y=log(median_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_50)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_50), y=log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_50), y=log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Haploid_threshold_100 <- c(0, checked_generations_vector_Haploid_threshold + 100)
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_Haploid_threshold_100)), span=1.5)
plot(x=(checked_generations_vector_Haploid_threshold_100), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Haploid_threshold_100), y=log(c(median_num_of_species_vector_Haploid_threshold) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_Haploid_threshold_100)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_Haploid_threshold_100), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Haploid_threshold_100), y=log(c(upper_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_Haploid_threshold_100)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_Haploid_threshold_100), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(c(lower_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("BDM", "Diploid threshold", "Haploid threshold"), col=c("black", "red", "orange"), fill=c("black", "red", "orange"))


## Haploid threshold, Haploid threshold assortative_mating

# log scale on y axis with range
# (plot Haploid log range.png, plot Haploid log range.eps)
# plot_haploid_log_range.pdf, plot_haploid_log_range.TIFF, plot_haploid_log_range.eps
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1), log(100+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50, 100)
lo <- loess((log(median_num_of_species_vector_Haploid_threshold + 1)) ~ (c(1, checked_generations_vector_Haploid_threshold)), span=1.5)
plot(x=c(1, checked_generations_vector_Haploid_threshold), y=predict(lo), xlab='generation', ylab='number of species', main='Haploid threshold & Haploid threshold assortative mating', cex.main=0.85, col='orange', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_Haploid_threshold, y=log(median_num_of_species_vector_Haploid_threshold + 1), xlab='generation', ylab='number of species', main='Haploid threshold & Haploid threshold assortative mating', cex.main=0.85, col='orange', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_Haploid_threshold + 1)) ~ (c(1, checked_generations_vector_Haploid_threshold)), span=1.5)
par(new=T)
plot(x=c(1, checked_generations_vector_Haploid_threshold), y=predict(lo), ann=F, col='orange', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Haploid_threshold, y=log(upper_num_of_species_vector_Haploid_threshold + 1), xlab='generation', ylab='number of species', main='Haploid threshold & Haploid threshold assortative mating', cex.main=0.85, col='orange', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Haploid_threshold + 1)) ~ (c(1, checked_generations_vector_Haploid_threshold)), span=1.5)
par(new=T)
plot(x=c(1, checked_generations_vector_Haploid_threshold), y=predict(lo), ann=F, col='orange', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_Haploid_threshold, y=log(lower_num_of_species_vector_Haploid_threshold + 1), xlab='generation', ylab='number of species', main='Haploid threshold & Haploid threshold assortative mating', cex.main=0.85, col='orange', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
par(new=T)

checked_generations_vector_Haploid_threshold_50 <- c(0, checked_generations_vector_Haploid_threshold + 50)
lo <- loess((log(median_num_of_species_vector_Haploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Haploid_threshold_50)), span=1.5)
plot(x=(checked_generations_vector_Haploid_threshold_50), y=predict(lo), ann=F, col='green', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_Haploid_threshold_50), y=log(median_num_of_species_vector_Haploid_threshold_assortative_mating + 1), ann=F, col='green', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Haploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Haploid_threshold_50)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_Haploid_threshold_50), y=predict(lo), ann=F, col='green', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Haploid_threshold_50), y=log(upper_num_of_species_vector_Haploid_threshold_assortative_mating + 1), ann=F, col='green', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
lo <- loess((log(lower_num_of_species_vector_Haploid_threshold_assortative_mating + 1)) ~ ((checked_generations_vector_Haploid_threshold_50)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_Haploid_threshold_50), y=predict(lo), ann=F, col='green', xlim=c(0, 4300), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_Haploid_threshold_50), y=log(lower_num_of_species_vector_Haploid_threshold_assortative_mating + 1), ann=F, col='green', xlim=c(0, 5000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="dotted")
legend("topleft", legend=c("Haploid threshold (5)", "Haploid threshold assortative mating (6)"), col=c("orange", "green"), fill=c("orange", "green"))


## BDM, BDM assortative_mating, Diploid threshold assortative_mating, Haploid threshold assortative_mating

# log scale on y axis with range
# (plot assortative mating and BDM log range.png, plot assortative mating and BDM log range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1), log(100+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50, 100)
lo <- loess((log(median_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold assortative mating \n& Haploid threshold assortative mating', cex.main=0.8, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(median_num_of_species_vector_BDM + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold assortative mating & Haploid threshold assortative mating', cex.main=0.75, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(upper_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(lower_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")

checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
par(new=T)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_50), y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")

checked_generations_vector_BDM_100 <- c(0, checked_generations_vector_BDM + 100)
par(new=T)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_100), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_150 <- c(0, (checked_generations_vector_BDM + 150)[1:5])
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_150)), span=1.5)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_150), y=log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_150)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_150)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
legend("topright", legend=c("BDM", "BDM assortative mating", "Diploid threshold assortative mating", "Haploid threshold assortative mating"), col=c("black", "blue", "purple", "green"), fill=c("black", "blue", "purple", "green"))


## All: BDM, BDM assortative_mating, Diploid threshold, Diploid threshold assortative_mating,
## Haploid threshold, Haploid threshold assortative_mating

# log scale on y axis with range
# (plot all log range.png, plot all log range.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating))))
exp(y_lim_max - 1)
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1), log(100+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50, 100)
lo <- loess((log(median_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), xlab='generation', ylab='number of species', main='All types', cex.main=0.85, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=checked_generations_vector_BDM, y=log(median_num_of_species_vector_BDM + 1), xlab='generation', ylab='number of species', main='BDM, BDM assortative mating, Diploid threshold assortative mating & Haploid threshold assortative mating', cex.main=0.75, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
lo <- loess((log(upper_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(upper_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_BDM + 1)) ~ (checked_generations_vector_BDM_0), span=0.4)
par(new=T)
plot(x=checked_generations_vector_BDM_0, y=predict(lo), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=checked_generations_vector_BDM, y=log(lower_num_of_species_vector_BDM + 1), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_50 <- c(0, checked_generations_vector_BDM + 50)
par(new=T)
lo <- loess((log(c(median_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_50), y=log(c(median_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(upper_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1)) ~ (checked_generations_vector_BDM_50), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_50), y=predict(lo), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM + 50), y=log(c(lower_num_of_species_vector_BDM_assortative_mating) + 1), ann=F, col='blue', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_100 <- c(0, checked_generations_vector_BDM + 100)
par(new=T)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_100), y=log(median_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(upper_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1)) ~ ((checked_generations_vector_BDM_100)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_100), y=predict(lo), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_100), y=log(lower_num_of_species_vector_Diploid_threshold[1:59] + 1), ann=F, col='red', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_150 <- c(0, checked_generations_vector_BDM + 150)
par(new=T)
lo <- loess((log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_150), y=log(median_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(upper_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1)) ~ ((checked_generations_vector_BDM_150)), span=0.4)
par(new=T)
plot(x=(checked_generations_vector_BDM_150), y=predict(lo), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_150), y=log(lower_num_of_species_vector_Diploid_threshold_assortative_mating[1:59] + 1), ann=F, col='purple', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_200 <- c(0, (checked_generations_vector_BDM + 200)[1:5])
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(median_num_of_species_vector_Haploid_threshold) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(upper_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold) + 1)) ~ ((checked_generations_vector_BDM_200)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_200), y=predict(lo), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_200), y=log(c(lower_num_of_species_vector_Haploid_threshold, rep(NA, (59-5))) + 1), ann=F, col='orange', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
par(new=T)

checked_generations_vector_BDM_250 <- c(0, (checked_generations_vector_BDM + 250)[1:5])
lo <- loess((log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
#plot(x=(checked_generations_vector_BDM_150), y=log(c(median_num_of_species_vector_Haploid_threshold_assortative_mating) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lwd=2)
lo <- loess((log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_250), y=log(c(upper_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
lo <- loess((log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating) + 1)) ~ ((checked_generations_vector_BDM_250)), span=1.5)
par(new=T)
plot(x=(checked_generations_vector_BDM_250), y=predict(lo), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=(checked_generations_vector_BDM_250), y=log(c(lower_num_of_species_vector_Haploid_threshold_assortative_mating, rep(NA, (59-5))) + 1), ann=F, col='green', xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", type="l", lty="longdash")
legend("topright", legend=c("BDM", "BDM assortative mating", "Diploid threshold", "Diploid threshold assortative mating", "Haploid threshold", "Haploid threshold assortative mating"), col=c("black", "blue", "red", "purple", "orange", "green"), fill=c("black", "blue", "red", "purple", "orange", "green"))





## Sympatric assortative mating

## (plot Sympatric assortative mating.png, plot Sympatric assortative mating.eps)
#y_lim_min <- 0.5
#(y_lim_max <- 1+max(c(upper_num_of_species_vector_sympatric_assortative_mating)))
#exp(y_lim_max - 1)
#y_tic_marks_location <- c(1, 2, 3)
#y_tic_marks_values <- c(1, 2, 3)

#plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, median_num_of_species_vector_sympatric_assortative_mating), xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), xlab='generation', ylab='number of species', main='Sympatric assortative mating', yaxt="n", xaxt="n", type="l", lwd=2)
#axis(1, at=c(0, 20000, 40000, 100000, 80000, 100000), labels=c("0", "20000", "40000", "100000", "80000", "100000"), las=1)                                     
#axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
#par(new=T)
#plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, upper_num_of_species_vector_sympatric_assortative_mating), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", xaxt="n", type="l", lty="longdash")
#par(new=T)
#plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, lower_num_of_species_vector_sympatric_assortative_mating), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", xaxt="n", type="l", lty="longdash")



# Paint by hand

# vector of checked generations
checked_generations_vector_sympatric_assortative_mating <- seq(999, 98999, 1000)
length(checked_generations_vector_sympatric_assortative_mating)

upper_num_of_species_vector_sympatric_assortative_mating <- c(
   1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1,
   2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1,
   1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2,
   1, 2, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 2,
   2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 3, 2, 1, 2, 1, 2, 1
)

median_num_of_species_vector_sympatric_assortative_mating <- c(
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
)

lower_num_of_species_vector_sympatric_assortative_mating <- c(
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
)

# (plot Sympatric assortative mating.png, plot Sympatric assortative mating.eps)
y_lim_min <- 0.5
(y_lim_max <- 1+max(c(upper_num_of_species_vector_sympatric_assortative_mating)))
exp(y_lim_max - 1)
y_tic_marks_location <- c(1, 2, 3)
y_tic_marks_values <- c(1, 2, 3)

plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, median_num_of_species_vector_sympatric_assortative_mating), xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), xlab='generation', ylab='number of species', main='Sympatric assortative mating', yaxt="n", xaxt="n", type="l", lwd=2, cex=2)
axis(1, at=c(0, 20000, 40000, 60000, 80000, 100000), labels=c("0", "20000", "40000", "60000", "80000", "100000"), las=1, cex=0.9)                                     
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2)
par(new=T)
plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, upper_num_of_species_vector_sympatric_assortative_mating), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", xaxt="n", type="l", lty="longdash")
par(new=T)
plot(x=c(0, checked_generations_vector_sympatric_assortative_mating), y=c(1, lower_num_of_species_vector_sympatric_assortative_mating), ann=F, xlim=c(0, 100000), ylim=c(y_lim_min, y_lim_max), yaxt="n", xaxt="n", type="l", lty="longdash")







######## boxplots of the number of genetic differences between two individuals
######## from different species in each of the cases (types of speciation)
######## for generations 5000 and 99000.


number_of_differences_vector_all <- c(number_of_differences_vector_BDM_5000_no_zeroes,
                                      number_of_differences_vector_BDM_assortative_mating_5000_no_zeroes,
                                      number_of_differences_vector_Diploid_threshold_assortative_mating_5000_no_zeroes,
                                      number_of_differences_vector_Haploid_threshold_5000_no_zeroes,
                                      number_of_differences_vector_Haploid_threshold_assortative_mating_5000_no_zeroes,
                                      
                                      number_of_differences_vector_BDM_99000_no_zeroes,
                                      number_of_differences_vector_BDM_assortative_mating_99000_no_zeroes,
                                      number_of_differences_vector_Diploid_threshold_assortative_mating_99000_no_zeroes)


case_vector_for_differences_all <- c(rep(1, length(number_of_differences_vector_BDM_5000_no_zeroes)),
                                     rep(2, length(number_of_differences_vector_BDM_assortative_mating_5000_no_zeroes)),
                                     rep(3, length(number_of_differences_vector_Diploid_threshold_assortative_mating_5000_no_zeroes)),
                                     rep(4, length(number_of_differences_vector_Haploid_threshold_5000_no_zeroes)),
                                     rep(5, length(number_of_differences_vector_Haploid_threshold_assortative_mating_5000_no_zeroes)),
                                     
                                     rep(8, length(number_of_differences_vector_BDM_99000_no_zeroes)),
                                     rep(9, length(number_of_differences_vector_BDM_assortative_mating_99000_no_zeroes)),
                                     rep(10, length(number_of_differences_vector_Diploid_threshold_assortative_mating_99000_no_zeroes)))


# (boxplot differences.png, boxplot differences.eps)
x_lim_min <- 0
x_lim_max <- 11
y_lim_min <- 0
y_lim_max <- 10000
x_tic_marks_location <- c(1:5, 8:10)
x_tic_marks_values <- c("BDM", "BDM \nassor. \nmating", "Dip. thres. \nassor. \nmating", "Hap. \nthres.", "Hap. thres. \nassor. \nmating", "BDM", "BDM \nassor. \nmating", "Dip. thres. \nassor. \nmating")
y_tic_marks_location <- c(0, 2000, 4000, 6000, 8000, 10000)
y_tic_marks_values <- c(0, 2000, 4000, 6000, 8000, 10000)
plot(x=NA, y=NA, xlab="generation 5000                                              generation 99000",
     ylab='number of genetic differences',
     main="Number of differences between individuals from different species",
     cex.main=0.85, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))
axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.5)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.5)

bplot(add=T, x=number_of_differences_vector_all, by=case_vector_for_differences_all,
      xlim=c(0,11),
      pos=c(1:5, 8:10),
      border=c("black", "blue", "purple", "orange", "green", "black", "blue", "purple"),
      ann=F, xaxt='n', yaxt='n',)


# log scale on y axis
# (boxplot differences log.png, boxplot differences log.eps)
x_lim_min <- 0
x_lim_max <- 11
y_lim_min <- 0.5
y_lim_max <- 1+log(10000)
x_tic_marks_location <- c(1:5, 8:10)
#x_tic_marks_values <- c("BDM (1)", "BDM \nassor. \nmating (2)", "Dip. thres. \nassor. \nmating (4)", "Hap. \nthres.(5)", "Hap. thres. \nassor. \nmating (6)", "BDM (1)", "BDM \nassor. \nmating (2)", "Dip. thres. \nassor. \nmating (4)")
x_tic_marks_values <- c("(1)", "(2)", "(4)", "(5)", "(6)", "(1)", "(2)", "(4)")
y_tic_marks_location <- c(log(1+1), log(2+1), log(5+1), log(10+1), log(20+1), log(50+1), log(100+1), log(200+1), log(500+1), log(1000+1), log(2000+1), log(5000+1), log(10000+1))
y_tic_marks_values <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
plot(x=NA, y=NA, xlab="generation 5000                                              generation 99000",
     ylab='number of genetic differences',
     main="Number of differences between individuals from different species",
     cex.main=1.2, xaxt='n', yaxt='n',
     xlim=c(x_lim_min, x_lim_max), ylim=c(y_lim_min, y_lim_max))
axis(1, at=x_tic_marks_location, labels=x_tic_marks_values, las=1, cex.axis=0.8)
axis(2, at=y_tic_marks_location, labels=y_tic_marks_values, las=2, cex.axis=0.75)
legend("topright", legend=c("BDM (1)", "BDM assortative mating (2)", "Diploid threshold assortative mating (4)", "Haploid threshold (5)", "Haploid threshold assortative mating (6)"), col=c("black", "blue", "purple", "orange", "green"), fill=c("black", "blue", "purple", "orange", "green"), cex = 0.85)



#bplot(add=T, x=log(number_of_differences_vector_all + 1), by=case_vector_for_differences_all,
#      xlim=c(0,11),
#      pos=c(1:5, 8:10),
#      border=c("black", "blue", "purple", "orange", "green", "black", "blue", "purple"),
#      ann=F, xaxt='n', yaxt='n',)

# Painting by hand

number_of_differences_vector_all <- c(167,  140, 167, 181, 181, 200,  9.5,  8, 13.5, 17, 26, 33,  
                                      3200, 9500, 9500, 9500, 9500, 9500, 9500, 9500, 9500, 9500,
                                      2500, 2650, 2650, 2800, 1580, 1940, 1940, 1940, 1940, 2075, 2075,  105, 105, 135, 175, 175 
                                      )
case_vector_for_differences_all <- c(1,  2, 2, 2, 2, 2,  3,  4, 4, 4, 4, 4,
                                     5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                                     6, 6, 6, 6,  7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8)

bplot(add=T, x=log(number_of_differences_vector_all), by=case_vector_for_differences_all,
      xlim=c(0,11),
      pos=c(1:5, 8:10),
      border=c("black", "blue", "purple", "orange", "green", "black", "blue", "purple"),
      ann=F, xaxt='n', yaxt='n',)






######## Are the species in BDM more stable than in BDM_assortative_mating?
######## That is, does an individual in the BDM case tend to stay in the same
######## species more than an individual in the BDM_assortative_mating case?
########
######## For each run, we count the total number of species changes for all
######## individuals.
######## We avoid counting the most rapid changes, which occur in the
######## BDM_assortative_mating case, by not counting generations with
######## over 8 different species. This makes the tests we perform more
######## conservative.



#### Count the total number of species changes for each run of the BDM case

# Matrix of first row of each generation.
# Each row is for a different generation.
# Each column is for a different run.
first_row_of_each_generation_inside_run__BDM <- matrix(rep(0, num_of_generations_checked__BDM * sample_size), ncol=sample_size)
for (current_run in 1:sample_size){
   current_gen <- 1
   current_row <- 1
   while (current_gen <= num_of_generations_checked__BDM){
      first_row_found <- FALSE
      while (!first_row_found){
         if (division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][current_row, "gen"] == current_gen){
            first_row_of_each_generation_inside_run__BDM[current_gen, current_run] <- current_row
            first_row_found <- TRUE
            current_gen <- current_gen + 1
         }
         current_row <- current_row + 1
      }
   }
}


# vector of total number of species changes in each generation of the BDM case
change_BDM <- rep(0, sample_size)


for (current_run in 1:sample_size){
   
   current_gen <- 2
   
   while (current_gen <= num_of_generations_checked__BDM){
   
      current_individual_to_check <- 1
      next_individual_to_check <- 1
      current_section_previous_generation <- 1
      current_section_current_generation <- 1
      
      backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
      backwards_found <- FALSE
      while (!backwards_found){
         if(is.na(division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ]
                  [first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run], "species number"])){
            backwards <- backwards + 1
         } else{
            backwards_found <- TRUE
         }
      }
      
      current_num_of_species_previous_generation <- 
         division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run], "species number"]
      current_num_of_species_current_generation <- 1
      last_individual_in_current_section_previous_generation <-
         division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run]), "section size"]
      last_individual_in_current_section_current_generation <-
         division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[current_gen, current_run]), "section size"]
      
      while (current_individual_to_check <= num_of_individuals_checked){
         next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                             as.numeric(last_individual_in_current_section_current_generation))
         if (division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run] + 
                                                                                                          current_section_previous_generation - 1), "species number"] !=
                division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen), current_run] + 
                                                                                                             current_section_current_generation - 1), "species number"]){
            change_BDM[current_run] <- change_BDM[current_run] + (next_individual_to_check - current_individual_to_check)
         }
         if (last_individual_in_current_section_previous_generation < 
                as.numeric(last_individual_in_current_section_current_generation)){
            current_num_of_species_previous_generation <- 
               division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run] + 
                                                                                                            current_section_previous_generation), "species number"]
            current_section_previous_generation <- current_section_previous_generation + 1
            last_individual_in_current_section_previous_generation <- 
               last_individual_in_current_section_previous_generation + 
               division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run] + 
                                                                                                            current_section_previous_generation - 1), "section size"]
         } else{
            if (last_individual_in_current_section_previous_generation >
                   as.numeric(last_individual_in_current_section_current_generation)){
               current_num_of_species_current_generation <- 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[current_gen, current_run] +
                                                                                                               current_section_current_generation), "species number"]
               current_section_current_generation <- current_section_current_generation + 1
               last_individual_in_current_section_current_generation <- 
                  last_individual_in_current_section_current_generation + 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[current_gen, current_run] + 
                                                                                                               current_section_current_generation - 1), "section size"]
            } else{
               current_num_of_species_previous_generation <- 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run] + 
                                                                                                               current_section_previous_generation), "species number"]
               current_section_previous_generation <- current_section_previous_generation + 1
               last_individual_in_current_section_previous_generation <- 
                  last_individual_in_current_section_previous_generation + 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[(current_gen - backwards), current_run] + 
                                                                                                               current_section_previous_generation - 1), "section size"]
               current_num_of_species_current_generation <- 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[current_gen, current_run] + 
                                                                                                               current_section_current_generation), "species number"]
               current_section_current_generation <- current_section_current_generation + 1
               last_individual_in_current_section_current_generation <- 
                  last_individual_in_current_section_current_generation + 
                  division_into_species_table_BDM[division_into_species_table_BDM[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM[current_gen, current_run] + 
                                                                                                               current_section_current_generation - 1), "section size"]
            }
         }
         current_individual_to_check <- next_individual_to_check
      } 
      
      current_gen <- current_gen + 1
      
   }
   
}


change_BDM
#  67  86  96  68   0 122  40 139 160  40



#### Count the total number of species changes for each run of the BDM_assortative_mating case

# Matrix of first row of each generation.
# Each row is for a different generation.
# Each column is for a different run.
first_row_of_each_generation_inside_run__BDM_assortative_mating <- matrix(rep(0, num_of_generations_checked__BDM_assortative_mating * sample_size), ncol=sample_size)
for (current_run in 1:sample_size){
   current_gen <- 1
   current_row <- 1
   while (current_gen <= num_of_generations_checked__BDM_assortative_mating){
      first_row_found <- FALSE
      while (!first_row_found){
         if (division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][current_row, "gen"] == current_gen){
            first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run] <- current_row
            first_row_found <- TRUE
            current_gen <- current_gen + 1
         }
         current_row <- current_row + 1
      }
   }
}


# vector of total number of species changes in each generation of the BDM_assortative_mating case
change_BDM_assortative_mating <- rep(0, sample_size)


for (current_run in 1:sample_size){
   
   current_gen <- 2
   
   while (current_gen <= num_of_generations_checked__BDM_assortative_mating){
         
      if (!is.na(division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ]
                 [first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run], "species number"])){
            
         current_individual_to_check <- 1
         next_individual_to_check <- 1
         current_section_previous_generation <- 1
         current_section_current_generation <- 1
         
         backwards <- 1  # How far backwards is the closest generation that doesn't have NA's for species numbers
         backwards_found <- FALSE
         while (!backwards_found){
            if(is.na(division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ]
                     [first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run], "species number"])){
               backwards <- backwards + 1
            } else{
               backwards_found <- TRUE
            }
         }
         
         current_num_of_species_previous_generation <- 
            division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run], "species number"]
         current_num_of_species_current_generation <- 1
         last_individual_in_current_section_previous_generation <-
            division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run]), "section size"]
         last_individual_in_current_section_current_generation <-
            division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run]), "section size"]
         
         while (current_individual_to_check <= num_of_individuals_checked){
            next_individual_to_check <- 1 + min(as.numeric(last_individual_in_current_section_previous_generation),
                                                as.numeric(last_individual_in_current_section_current_generation))
            if (division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run] + 
                                                                                                                                                   current_section_previous_generation - 1), "species number"] !=
                   division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen), current_run] + 
                                                                                                                                                      current_section_current_generation - 1), "species number"]){
               change_BDM_assortative_mating[current_run] <- change_BDM_assortative_mating[current_run] + (next_individual_to_check - current_individual_to_check)
            }
            if (last_individual_in_current_section_previous_generation < 
                   as.numeric(last_individual_in_current_section_current_generation)){
               current_num_of_species_previous_generation <- 
                  division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run] + 
                                                                                                                                                     current_section_previous_generation), "species number"]
               current_section_previous_generation <- current_section_previous_generation + 1
               last_individual_in_current_section_previous_generation <- 
                  last_individual_in_current_section_previous_generation + 
                  division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run] + 
                                                                                                                                                     current_section_previous_generation - 1), "section size"]
            } else{
               if (last_individual_in_current_section_previous_generation >
                      as.numeric(last_individual_in_current_section_current_generation)){
                  current_num_of_species_current_generation <- 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run] +
                                                                                                                                                        current_section_current_generation), "species number"]
                  current_section_current_generation <- current_section_current_generation + 1
                  last_individual_in_current_section_current_generation <- 
                     last_individual_in_current_section_current_generation + 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run] + 
                                                                                                                                                        current_section_current_generation - 1), "section size"]
               } else{
                  current_num_of_species_previous_generation <- 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run] + 
                                                                                                                                                        current_section_previous_generation), "species number"]
                  current_section_previous_generation <- current_section_previous_generation + 1
                  last_individual_in_current_section_previous_generation <- 
                     last_individual_in_current_section_previous_generation + 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[(current_gen - backwards), current_run] + 
                                                                                                                                                        current_section_previous_generation - 1), "section size"]
                  current_num_of_species_current_generation <- 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run] + 
                                                                                                                                                        current_section_current_generation), "species number"]
                  current_section_current_generation <- current_section_current_generation + 1
                  last_individual_in_current_section_current_generation <- 
                     last_individual_in_current_section_current_generation + 
                     division_into_species_table_BDM_assortative_mating[division_into_species_table_BDM_assortative_mating[, "run"]==current_run, ][(first_row_of_each_generation_inside_run__BDM_assortative_mating[current_gen, current_run] + 
                                                                                                                                                        current_section_current_generation - 1), "section size"]
               }
            }
            current_individual_to_check <- next_individual_to_check
         } 
         
      }
      
      current_gen <- current_gen + 1
         
   }
   
}


change_BDM_assortative_mating
#  777 865 531 632 857 442 793 874 695 821



#### Test for normallity using the Shapiro-Wilk normality test.

shapiro.test(change_BDM)
#  W = 0.98, p-value = 0.965   =>   Can be considered normal.

shapiro.test(change_BDM_assortative_mating)
#  W = 0.8779, p-value = 0.1235   =>   Can be considered normal.



#### Levene's test for homoscedasticity

change_BDM_BDM_assortative_mating <- c(change_BDM, change_BDM_assortative_mating)
group_change_BDM_BDM_assortative_mating <- as.factor(c(rep(1, length(change_BDM)), rep(2, length(change_BDM_assortative_mating))))
leveneTest(y=change_BDM_BDM_assortative_mating, group=group_change_BDM_BDM_assortative_mating)
#  F value = 4.4117, Pr(>F) = 0.05005   =>  
#  Can be considered homoscedastic, but borderline.



#### Test the alternative hypothesis that species are more stable in the BDM case
#### (i.e., the number of total changes in species is lower) than in the 
#### BDM_assortative_mating case.
#### To be conservative, make it a two-sided test.
#### Use the non-parametric two sample Kruskal?Wallis two-sided test for the
#### null hypothesis that the two samples originate from the same distribution,
#### against the alternative that one sample stochastically dominates the other.


median(change_BDM)
#  77

sqrt(var(change_BDM))
#  49.23594


median(change_BDM_assortative_mating)
#  785

sqrt(var(change_BDM_assortative_mating))
#  150.1503


kruskal.test(list(change_BDM, change_BDM_assortative_mating))
# chi-squared = 14.2965, df = 1, p-value = 0.0001562   =>   
# The number of total changes in species in BDM < BDM_assortative_mating.
