library(tidyverse)

Power.Overlap.values <- read.csv("Vals.csv")[,2:11]

ratio <- 2/c(8, 4, 2, 1, 1/2)

threshold.vec <- Power.Overlap.values[1,]

Power <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[2:6,])))
Overlap <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[7:11,])))

colnames(Power) <- c("threshold", "ratio1", "ratio2",
                     "ratio3", "ratio4", "ratio5")

colnames(Overlap) <- c("threshold", "ratio1", "ratio2",
                    "ratio3", "ratio4", "ratio5")
power.plot <- Power %>%
  ggplot(aes(x = threshold, y = ratio1)) + 
  geom_vline(xintercept = 1, linetype = "dotted", size = 0.5) +
  geom_line(aes(y = ratio1, color = "ratio1"), linetype = "solid") + 
  geom_point(aes(y = ratio1, color = "ratio1")) +
  geom_line(aes(y = ratio2, color = "ratio2"), linetype = "solid") + 
  geom_point(aes(y = ratio2, color = "ratio2")) +
  geom_line(aes(y = ratio3, color = "ratio3"), linetype = "solid") + 
  geom_point(aes(y = ratio3, color = "ratio3")) +
  geom_line(aes(y = ratio4, color = "ratio4"), linetype = "solid") + 
  geom_point(aes(y = ratio4, color = "ratio4")) +
  geom_line(aes(y = ratio5, color = "ratio5"), linetype = "solid") +
  geom_point(aes(y = ratio5, color = "ratio5")) +
  scale_color_manual(breaks = c("ratio1", "ratio2",
                                "ratio3", "ratio4", "ratio5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = ratio) +
  labs(title = "Empirical Power", x = "Combined Signal Strength",
       y = "Power", color = "Ratio") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85,0.5),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = power.plot, filename = "ratio_power.pdf", device = "pdf",
       width = 4, height = 3)

overlap.plot <- Overlap %>%
  ggplot(aes(x = threshold, y = ratio1)) + 
  geom_vline(xintercept = 1, linetype = "dotted", size = 0.5) +
  geom_line(aes(y = ratio1, color = "ratio1"), linetype = "solid") + 
  geom_point(aes(y = ratio1, color = "ratio1")) +
  geom_line(aes(y = ratio2, color = "ratio2"), linetype = "solid") + 
  geom_point(aes(y = ratio2, color = "ratio2")) +
  geom_line(aes(y = ratio3, color = "ratio3"), linetype = "solid") + 
  geom_point(aes(y = ratio3, color = "ratio3")) +
  geom_line(aes(y = ratio4, color = "ratio4"), linetype = "solid") + 
  geom_point(aes(y = ratio4, color = "ratio4")) +
  geom_line(aes(y = ratio5, color = "ratio5"), linetype = "solid") +
  geom_point(aes(y = ratio5, color = "ratio5")) +
  scale_color_manual(breaks = c("ratio1", "ratio2",
                                "ratio3", "ratio4", "ratio5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = ratio) +
  labs(title = "Empirical Overlap", x = "Combined Signal Strength",
       y = "Overlap", color = "Ratio") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.15,0.6),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = overlap.plot, filename = "ratio_overlap.pdf", device = "pdf",
       width = 4, height = 3)

