1180 - 281 / 282 / 288

2638 - 905 / 924 / 861

3818 - 1360 / 1180 / 1292

4642 - 1752 / 1599 / 1458

5276 - 2168 / 2146 / 2086

# MOUSE

9190 - 200 / 205 / 215

13630 - 387 / 395 / 401

16468 - 516 / 528 / 519

28507 - 1304 / 1238 / 1377

37402 - 2177 / 2260 / 2171

44879 - 2864 / 2701 / 3004

50296 - 4218 / 4962 / 4195

61637 - 5733 / 5948 / 5350

library(ggplot2)

noCell <- rep(c(1180, 2638,3818, 4642, 5276,
                9190,13630,16468,28507,37402,44879,50296,61637), each = 3)

times <- c(
  281, 282, 288,
  905,924,861,
  1360,1180,1292,
  1752,1599,1458,
  2168,2146,2086,

  200,205,215,
  387,395,401,
  516,528,519,
  1304,1238,1377,
  2177,2260,2171,
  2864,2701,3004,
  4218,4962,4195,
  5733,5948,5350
)
species <- c(rep('Human - Curated', 3*5), rep('Mouse - KEGG', 3*8))
tt <- data.frame(noCell, times, species )

ggplot(tt, aes(x = noCell, y = times, shape = species, color = species)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme(
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position =  'none'
  ) +
  facet_wrap(~ species, scales = 'free')+
  ylab('Running Time ( Seconds )') +
  xlab('Number of Cells')



