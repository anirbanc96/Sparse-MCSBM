library(tidyverse)
require(latex2exp)

Power.Overlap.values <- as_tibble(read_csv("Vals.csv"))

corr <- c(0, 0.001, 0.01, 0.1)

threshold.vec <- Power.Overlap.values[1,]

Power <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[2:5,])))
Overlap <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[6:9,])))

colnames(Power) <- c("threshold", "corr1", "corr2",
                     "corr3", "corr4")

colnames(Overlap) <- c("threshold", "corr1", "corr2",
                       "corr3", "corr4")
Power %>%
  ggplot(aes(x = threshold, y = corr1)) + 
  geom_line(aes(y = corr1, color = "corr1"), linetype = "dashed") + 
  geom_point(aes(y = corr1, color = "corr1")) +
  geom_line(aes(y = corr2, color = "corr2"), linetype = "dashed") + 
  geom_point(aes(y = corr2, color = "corr2")) +
  geom_line(aes(y = corr3, color = "corr3"), linetype = "dashed") + 
  geom_point(aes(y = corr3, color = "corr3")) +
  geom_line(aes(y = corr4, color = "corr4"), linetype = "dashed") + 
  geom_point(aes(y = corr4, color = "corr4")) +
  scale_color_manual(breaks = c("corr1", "corr2",
                                "corr3", "corr4"),
                     values = c("cornflowerblue", "chocolate1", "deeppink",
                                "forestgreen"),
                     labels = c(expression(paste(rho, " = ", 0)),
                                expression(paste(rho, " = ", 0.001)),
                                expression(paste(rho, " = ", 0.01)),
                                expression(paste(rho, " = ", 0.1)))) +
  labs(title = "Power of Belief Propagation", x = "Threshold",
       y = "Power", color = "Correlation") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.75,0.5),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = power.plot, filename = "corr_power.pdf", device = "pdf",
       width = 4, height = 3)

overlap.plot <- Overlap %>%
  ggplot(aes(x = threshold, y = corr1)) + 
  geom_line(aes(y = corr1, color = "corr1"), linetype = "dashed") + 
  geom_point(aes(y = corr1, color = "corr1")) +
  geom_line(aes(y = corr2, color = "corr2"), linetype = "dashed") + 
  geom_point(aes(y = corr2, color = "corr2")) +
  geom_line(aes(y = corr3, color = "corr3"), linetype = "dashed") + 
  geom_point(aes(y = corr3, color = "corr3")) +
  geom_line(aes(y = corr4, color = "corr4"), linetype = "dashed") + 
  geom_point(aes(y = corr4, color = "corr4")) +
  scale_color_manual(breaks = c("corr1", "corr2",
                                "corr3", "corr4"),
                     values = c("cornflowerblue", "chocolate1", "deeppink",
                                "forestgreen"),
                     labels = c(expression(paste(rho, " = ", 0)),
                                expression(paste(rho, " = ", 0.001)),
                                expression(paste(rho, " = ", 0.01)),
                                expression(paste(rho, " = ", 0.1)))) +
  labs(title = "Overlap of Belief Propagation", x = "Threshold",
       y = "Power", color = "Correlation") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.2,0.5),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = overlap.plot, filename = "corr_overlap.pdf", device = "pdf",
       width = 4, height = 3)

