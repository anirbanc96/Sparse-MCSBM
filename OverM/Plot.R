library(tidyverse)

Power.Overlap.values <- as_tibble(read_csv("Vals.csv"))

m <- c(1, 3, 5, 7, 10)

threshold.vec <- Power.Overlap.values[1,]

Power <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[2:6,])))
Overlap <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[7:11,])))

colnames(Power) <- c("threshold", "m1", "m2",
                     "m3", "m4", "m5")

colnames(Overlap) <- c("threshold", "m1", "m2",
                       "m3", "m4", "m5")
power.plot <- Power %>%
  ggplot(aes(x = threshold, y = m1)) + 
  geom_line(aes(y = m1, color = "m1"), linetype = "dashed") + 
  geom_point(aes(y = m1, color = "m1")) +
  geom_line(aes(y = m2, color = "m2"), linetype = "dashed") + 
  geom_point(aes(y = m2, color = "m2")) +
  geom_line(aes(y = m3, color = "m3"), linetype = "dashed") + 
  geom_point(aes(y = m3, color = "m3")) +
  geom_line(aes(y = m4, color = "m4"), linetype = "dashed") + 
  geom_point(aes(y = m4, color = "m4")) +
  geom_line(aes(y = m5, color = "m5"), linetype = "dashed") +
  geom_point(aes(y = m5, color = "m5")) +
  scale_color_manual(breaks = c("m1", "m2",
                                "m3", "m4", "m5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = m) +
  labs(title = "Power of Belief Propagation", x = "Threshold",
       y = "Power", color = "m") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85,0.5),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = power.plot, filename = "m_power.pdf", device = "pdf",
       width = 4, height = 3)

overlap.plot <- Overlap %>%
  ggplot(aes(x = threshold, y = m1)) + 
  geom_line(aes(y = m1, color = "m1"), linetype = "dashed") + 
  geom_point(aes(y = m1, color = "m1")) +
  geom_line(aes(y = m2, color = "m2"), linetype = "dashed") + 
  geom_point(aes(y = m2, color = "m2")) +
  geom_line(aes(y = m3, color = "m3"), linetype = "dashed") + 
  geom_point(aes(y = m3, color = "m3")) +
  geom_line(aes(y = m4, color = "m4"), linetype = "dashed") + 
  geom_point(aes(y = m4, color = "m4")) +
  geom_line(aes(y = m5, color = "m5"), linetype = "dashed") +
  geom_point(aes(y = m5, color = "m5")) +
  scale_color_manual(breaks = c("m1", "m2",
                                "m3", "m4", "m5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = m) +
  labs(title = "Overlap from Belief Propagation", x = "Threshold",
       y = "Overlap", color = "m") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.15,0.6),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = overlap.plot, filename = "m_overlap.pdf", device = "pdf",
       width = 4, height = 3)

