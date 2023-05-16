rm(list=ls())
source("helper_functions.R")
source("simulation_functions.R")


# caterpillar

n_rep = 100
n = 16
h = 500
seq_length <- seq(100, 700, 100)

set.seed(1)
cater <- Get_result_binary_tree_gaussian(n_rep, n, h, seq_length, type = 'caterpillar')

mean <- as.vector(t(matrix(colMeans(cater), ncol = length(seq_length))))
sd <- as.vector(t(matrix(apply(cater, 2, sd), ncol = length(seq_length))))
method <- rep(c("NJ", "ECR"), each = length(seq_length))
sequence_len <- rep(seq_length, 2)
cater_df <- data.frame(mean = mean, sd = sd, sequence_len = sequence_len, method = method)
cater_df$tree <- rep("caterpillar", nrow(cater_df))

# ggplot(cater_df, aes(x=sequence_len, y=mean, colour=method, shape = method)) + 
#   geom_line() +
#   geom_point() +
#   labs(y= "RF distance", x = "number of samples used to estimate D")+
#   scale_color_manual(values=c('#999999','#E69F00'))+
#   labs(y= "RF distance", x = "sequence length")+
#   theme_classic()


# balanced binary

n_rep = 100
n = 16
h = 500
seq_length <- seq(100, 700, 100)

set.seed(1)
balance <- Get_result_binary_tree_gaussian(n_rep, n, h, seq_length, type = 'balanced')

mean <- as.vector(t(matrix(colMeans(balance), ncol = length(seq_length))))
sd <- as.vector(t(matrix(apply(balance, 2, sd), ncol = length(seq_length))))
method <- rep(c("NJ", "ECR"), each = length(seq_length))
sequence_len <- rep(seq_length, 2)
balance_df <- data.frame(mean = mean, sd = sd, sequence_len = sequence_len, method = method)
balance_df$tree <- rep("balanced binary", nrow(balance_df))

# ggplot(balance_df, aes(x=sequence_len, y=mean, colour=method, shape = method)) + 
#   geom_line() +
#   geom_point() +
#   labs(y= "RF distance", x = "number of samples used to estimate D")+
#   scale_color_manual(values=c('#999999','#E69F00'))+
#   labs(y= "RF distance", x = "sequence length")+
#   theme_classic()

# random non-binary

n_rep = 100
n = 16
h = 500
seq_length <- round(exp(seq(4, 10)))

set.seed(1)
random_non_binary <- Get_result_non_binary_tree_gaussian(n_rep, n, h, seq_length)

mean <- as.vector(t(matrix(colMeans(random_non_binary), nrow = 4)))
sd <- as.vector(t(matrix(apply(random_non_binary, 2, sd), nrow = 4)))
method <- rep(c("NJ", "NJ", "ECR", "ECR"), each = 7)
sequence_len <- rep(seq_length, 4)
random_non_binary_df <- data.frame(mean = mean, sd = sd, sequence_len = sequence_len, method = method)
random_non_binary_df$tree <- rep(c("non-binary", "non-binary (k unknown)", "non-binary", "non-binary (k unknown)"), each = 7)

# ggplot(random_non_binary_df, aes(x=sequence_len, y=mean, colour=method, shape = method)) + 
#   geom_line() +
#   geom_point() +
#   facet_wrap(~tree, scales = "free_x")+ 
#   labs(y= "RF distance", x = "sequence length")+
#   scale_x_continuous(trans='log2') +
#   scale_color_manual(values=c('#999999','#E69F00'))+
#   theme_classic()


# combine into one dataframe
all <- rbind(cater_df, balance_df, random_non_binary_df)
saveRDS(all, file = "data/all_simu_results.RDS")


# plot

library(ggplot2)
library(ggh4x)

pdf(file = "plots/binary_cater_multi_gaussian.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 2)
ggplot(all, aes(x=sequence_len, y=mean, colour=method, shape = method)) + 
  geom_line() +
  geom_point() +
  facet_wrap(~tree, scales = "free_x", nrow = 1)+
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position=position_dodge(0.5))+
  facetted_pos_scales(
    x = list(
      tree == "non-binary" ~ scale_x_continuous(trans='log2'),
      tree == "non-binary (k unknown)" ~ scale_x_continuous(trans='log2')
    )
  )+
  labs(y= "RF distance", x = "number of samples used to estimate D")+
  #scale_x_continuous(trans='log2') +
  scale_color_manual(values=c('#999999','#E69F00'))+
  theme_classic()+
  theme( axis.text = element_text( size = 10 ),
         axis.text.x = element_text( size = 10 ),
         axis.title = element_text( size = 10 ),
         legend.title = element_text(size=8),
         legend.text = element_text(size=8),
         strip.text = element_text(size = 9, face = "bold"))
dev.off()

# error bar 

t(spread(all[all$tree=='caterpillar',-c(1, 5)], method, sd))

t(spread(all[all$tree=='balanced binary',-c(1, 5)], method, sd))

t(spread(all[all$tree=='non-binary',-c(1, 5)], method, sd))

t(spread(all[all$tree=='non-binary (k unknown)',-c(1, 5)], method, sd))







