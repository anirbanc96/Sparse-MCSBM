library(tidyverse)
require(latex2exp)

Power.Overlap.values <- as_tibble(read_csv("Vals.csv"))

n <- c(200, 400, 600, 800, 1000)

p <- c(250, 500, 750, 1000, 1250)

threshold.vec <- Power.Overlap.values[1,]

Power <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[2:6,])))
Overlap <- as.data.frame(t(rbind(threshold.vec, Power.Overlap.values[7:11,])))

colnames(Power) <- c("threshold", "n1", "n2",
                     "n3", "n4", "n5")

colnames(Overlap) <- c("threshold", "n1", "n2",
                       "n3", "n4", "n5")
power.plot <- Power %>%
  ggplot(aes(x = threshold, y = n1)) + 
  geom_vline(xintercept = 1, linetype = "dotted", size = 0.5) +
  geom_line(aes(y = n1, color = "n1"), linetype = "solid") + 
  geom_point(aes(y = n1, color = "n1")) +
  geom_line(aes(y = n2, color = "n2"), linetype = "solid") + 
  geom_point(aes(y = n2, color = "n2")) +
  geom_line(aes(y = n3, color = "n3"), linetype = "solid") + 
  geom_point(aes(y = n3, color = "n3")) +
  geom_line(aes(y = n4, color = "n4"), linetype = "solid") + 
  geom_point(aes(y = n4, color = "n4")) +
  geom_line(aes(y = n5, color = "n5"), linetype = "solid") +
  geom_point(aes(y = n5, color = "n5")) +
  scale_color_manual(breaks = c("n1", "n2",
                                "n3", "n4", "n5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = c(paste(TeX(r'(n = )'), n[1], ",", TeX(r'(p = )'), p[1]), 
                                paste(TeX(r'(n = )'), n[2], ",", TeX(r'(p = )'), p[2]),
                                paste(TeX(r'(n = )'), n[3], ",", TeX(r'(p = )'), p[3]),
                                paste(TeX(r'(n = )'), n[4], ",", TeX(r'(p = )'), p[4]),
                                paste(TeX(r'(n = )'), n[5], ",", TeX(r'(p = )'), p[5]))) +
  labs(title = "Empirical Power", x = "Combined Signal Strength",
       y = "Power", color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.75,0.5),
        legend.title=element_blank(),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = power.plot, filename = "n_power.pdf", device = "pdf",
       width = 4, height = 3)

overlap.plot <- Overlap %>%
  ggplot(aes(x = threshold, y = n1)) + 
  geom_vline(xintercept = 1, linetype = "dotted", size = 0.5) +
  geom_line(aes(y = n1, color = "n1"), linetype = "solid") + 
  geom_point(aes(y = n1, color = "n1")) +
  geom_line(aes(y = n2, color = "n2"), linetype = "solid") + 
  geom_point(aes(y = n2, color = "n2")) +
  geom_line(aes(y = n3, color = "n3"), linetype = "solid") + 
  geom_point(aes(y = n3, color = "n3")) +
  geom_line(aes(y = n4, color = "n4"), linetype = "solid") + 
  geom_point(aes(y = n4, color = "n4")) +
  geom_line(aes(y = n5, color = "n5"), linetype = "solid") +
  geom_point(aes(y = n5, color = "n5")) +
  scale_color_manual(breaks = c("n1", "n2",
                                "n3", "n4", "n5"),
                     values = c("cornflowerblue", "darkred", "deeppink",
                                "forestgreen","chocolate1"),
                     labels = c(paste(TeX(r'(n = )'), n[1], ",", TeX(r'(p = )'), p[1]), 
                                paste(TeX(r'(n = )'), n[2], ",", TeX(r'(p = )'), p[2]),
                                paste(TeX(r'(n = )'), n[3], ",", TeX(r'(p = )'), p[3]),
                                paste(TeX(r'(n = )'), n[4], ",", TeX(r'(p = )'), p[4]),
                                paste(TeX(r'(n = )'), n[5], ",", TeX(r'(p = )'), p[5]))) +
  labs(title = "Empirical Overlap", x = "Combined Signal Strength",
       y = "Average Overlap", color = "Legend Title\n") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.75,0.32),
        legend.margin=margin(c(0,3,2,2)),
        legend.title=element_blank(),
        legend.text.align = 0,
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = overlap.plot, filename = "n_overlap.pdf", device = "pdf",
       width = 4, height = 3)

