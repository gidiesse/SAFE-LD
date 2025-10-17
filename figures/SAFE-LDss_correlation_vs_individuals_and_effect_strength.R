## Script to create figure 2.
## This one is annoying because since experiments at 1%, 5%, 10% and 20%  variance
## were distinct, I simply copied the values of correlation achieved manually
## to group them by number of individuals. There is no more efficient way
## I could think of doing this. 
## So inputs are:
## number of indivudals, % of variance explained and correlation value between
## LD and SAFE-LD

library(ggplot2)

#### NON-NULL graph 200 individuals 

# Sample data: 2 methods, 5 x-values each
df_200 <- data.frame(
  x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
  y = c(0.9931585, 0.9931551, 0.9931611, 0.9933897, 0.9931477,
        0.9933101, 0.9931215, 0.9935193, 0.9932163, 0.9932618,
        0.9932038, 0.9931508, 0.9931271, 0.9930174, 0.9930088, 
        0.9932453, 0.9931546, 0.9929744, 0.9927907, 0.9926730),
  sd = sqrt(c(2.408760e-07, 2.517760e-07, 1.112390e-07, 1.683756e-07, 1.868957e-07,
              1.183049e-07, 1.203240e-07, 1.855860e-07, 2.318222e-07, 7.822464e-08,
              2.972457e-07, 3.750702e-07, 1.284857e-07, 5.967413e-07, 2.653067e-07, 
              1.278406e-07, 5.213642e-07, 1.869197e-07, 2.733850e-07, 3.747591e-07)),
  method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
  Individuals = 200
)

# Convert x to factor if needed for clean spacing
df_200$x <- factor(df_200$x)

# Plot
ggplot(df_200, aes(x = x, y = y, color = method, group = method)) +
  #geom_line() +
  geom_point(position=position_dodge(width = .5, preserve = "total"),size=3) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.2, size = 0.3,position=position_dodge(width = .5, preserve = "total"),size=5) +
  #geom_line(linewidth = 0.8) +
  coord_cartesian(ylim = c(0.99, 1))  +
  labs(x = "Percentage of phenotypes with effect", y = "Average correlation across 10 iterations", title = "Comparison of Different effect size for 200 individuals - 10 iterations - threshold 1e-6") +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#4285f4","#34a853","#fbbc05","#ea4335")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


#### NON-NULL graph 400 individuals 

# Sample data: 2 methods, 5 x-values each
df_400 <- data.frame(
  x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
  y = c(0.9926162, 0.9924265, 0.9927486, 0.9924438, 0.9926797, 
        0.9927022, 0.9924535, 0.9928974, 0.9926998, 0.9922994, 
        0.9926093, 0.9924849, 0.9927372, 0.9923418, 0.9916486, 
        0.9925326, 0.9926620, 0.9922661, 0.9921320, 0.9919663),
  sd = sqrt(c(1.745802e-07, 1.952433e-07, 7.855357e-08, 1.649013e-07, 1.348411e-07,
              1.710704e-07, 1.609897e-07, 2.077961e-07, 1.771022e-07, 2.626636e-07, 
              1.871847e-07, 2.552364e-07, 2.294663e-07, 2.198230e-07, 1.723493e-07, 
              1.032431e-07, 1.713187e-07, 1.517595e-07, 1.791233e-07, 3.224466e-07)),
  method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
  Individuals = 400
)

# Convert x to factor if needed for clean spacing
df_400$x <- factor(df_400$x)

# Plot
ggplot(df_400, aes(x = x, y = y, color = method, group = method)) +
  #geom_line() +
  geom_point(position=position_dodge(width = .5, preserve = "total"),size=3) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.2, size = 0.3,position=position_dodge(width = .5, preserve = "total"),size=5) +
  #geom_line(linewidth = 0.8) +
  coord_cartesian(ylim = c(0.99, 1))  +
  labs(x = "Percentage of phenotypes with effect", y = "Average correlation across 10 iterations", title = "Comparison of Different effect size for 400 individuals - 10 iterations - threshold 1e-6") +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#4285f4","#34a853","#fbbc05","#ea4335")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
#### NON-NULL graph 1000 individuals 

# Sample data: 2 methods, 5 x-values each
df_1000 <- data.frame(
  x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
  y = c(0.9923370, 0.9923325, 0.9924902, 0.9922653, 0.9921635,
        0.9922110, 0.9923147, 0.9920592, 0.9918050, 0.9912965, 
        0.9920040, 0.9923185, 0.9922831, 0.9918394, 0.9913140,
        0.9922053, 0.9922453, 0.9921136, 0.9918589, 0.9915898),
  sd = sqrt(c(6.955446e-08, 1.553720e-07, 1.349810e-07, 2.698378e-07, 2.149716e-07,
              2.891483e-07, 7.398617e-08, 1.940100e-07, 2.965363e-07, 2.560318e-07,
              1.832391e-07, 3.589883e-07, 2.104146e-07, 2.569266e-07, 6.909829e-07,
              3.598396e-07, 2.719519e-07, 2.121434e-07, 2.551806e-07, 3.056797e-07
              
              )),
  method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
  Individuals=1000
)

# Convert x to factor if needed for clean spacing
df_1000$x <- factor(df_1000$x)

# Plot
ggplot(df_1000, aes(x = x, y = y, color = method, group = method)) +
  #geom_line() +
  geom_point(position=position_dodge(width = .5, preserve = "total"),size=3) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.2, size = 0.3,position=position_dodge(width = .5, preserve = "total"),size=5) +
  #geom_line(linewidth = 0.8) +
  coord_cartesian(ylim = c(0.99, 1))  +
  labs(x = "Percentage of phenotypes with effect", y = "Average correlation across 10 iterations", title = "Comparison of Different effect size for 1000 individuals - 10 iterations - threshold 1e-6") +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = c("#4285f4","#34a853","#fbbc05","#ea4335")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
# what if I put them all together?
df_all <- rbind(df_200, df_400, df_1000)


#####
# number of individuals phenotypes (e.g., 200, 500, 1000, 5000…) on x-axis and average correlation on y-axis and effect 
# (e.g., 0.1%, 0.5%, 1%, 5%, 10%)  with colored lines, and % of phenotypes with effect with facets (e.g., 0.1%, 0.5%, 1%, 5%, 10%)

# df_0.1_perc <- data.frame(
#   x = rep(c(200, 400, 1000), times = 5),
#   y = c(0.9931585, 0.9926162, 0.9923370, 
#         0.9933101, 0.9927022, 0.9922110, 
#         0.9932038, 0.9926093, 
#         0.9932453, 0.9925326, ),
#   sd = sqrt(c(2.408760e-07, 1.745802e-07, 6.955446e-08, 
#               1.183049e-07, 1.710704e-07, 
#               2.972457e-07, 1.871847e-07, 
#               1.278406e-07, 1.032431e-07, )),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   percentage_phenos_with_effect = 0.1
# )
# 
# df_0.5_perc <- data.frame(
#   x = rep(c(200, 400, 1000), times = 5),
#   y = c(0.9931551, 0.9924265, 0.9923325,
#         0.9931215, 0.9924535, 
#         0.9931508, 0.9924849, 
#         0.9931546, 0.9926620, 
#          ),
#   sd = sqrt(c(2.517760e-07, 1.952433e-07, 1.553720e-07, 
#               1.203240e-07, 1.609897e-07, 
#               3.750702e-07, 2.552364e-07, 
#               5.213642e-07, 1.713187e-07, 
#                )),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   percentage_phenos_with_effect = 0.1
# )
# 
# df_1_perc <- data.frame(
#   x = rep(c(200, 400, 1000), times = 5),
#   y = c(0.9931611, 0.9927486, 0.9924902, 
#         0.9935193, 0.9928974, 
#         0.9931271, 0.9927372, 
#         0.9929744, 0.9922661, 
#          ),
#   sd = sqrt(c(1.112390e-07, 7.855357e-08, 1.349810e-07, 
#               1.855860e-07, 2.077961e-07, 
#               1.284857e-07, 2.294663e-07, 
#               1.869197e-07, 1.517595e-07, 
#                )),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   percentage_phenos_with_effect = 0.1
# )
# 
# df_5_perc <- data.frame(
#   x = rep(c(200, 400, 1000), times = 5),
#   y = c(0.9933897, 0.9924438, 0.9922653, 
#         0.9932163, 0.9926998, 
#         0.9930174, 0.9923418, 
#         0.9927907, 0.9921320, ),
#   sd = sqrt(c(1.683756e-07, 1.649013e-07, 2.698378e-07, 
#               2.318222e-07, 1.771022e-07, 
#               5.967413e-07, 2.198230e-07, 
#               2.733850e-07, 1.791233e-07, )),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   percentage_phenos_with_effect = 0.1
# )
# 
# df_10_perc <- data.frame(
#   x = rep(c(200, 400, 1000), times = 5),
#   y = c(0.9931477, 0.9926797, 0.9921635, 
#         0.9932618, 0.9922994, 
#         0.9930088, 0.9916486, 
#         0.9926730, 0.9919663, ),
#   sd = sqrt(c(1.868957e-07, 1.348411e-07, 2.149716e-07, 
#               7.822464e-08, 2.626636e-07, 
#               2.653067e-07, 1.723493e-07, 
#               3.747591e-07, 3.224466e-07, )),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   percentage_phenos_with_effect = 0.1
# )
# 
# #######
# 
# df_200_8 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9931717, 0.9931746, 0.9931732, 0.9933891, 0.9931812,
#         0.9933236, 0.9931318, 0.9935545, 0.9932468, 0.9932982,
#         0.9932204, 0.9931928, 0.9931817, 0.9930886,0.9931776, 
#         0.9932534, 0.9931895, 0.9929726, 0.9928621, 0.9927849 ),
#   sd = sqrt(c(2.367412e-07, 2.459882e-07, 1.074160e-07, 1.662042e-07, 1.616948e-07,
#               1.161549e-07, 1.156980e-07, 1.913348e-07, 2.125889e-07, 7.136778e-08,
#               2.806832e-07, 3.489014e-07, 1.459871e-07, 5.015743e-07, 2.378483e-07, 
#               1.409116e-07, 5.072734e-07, 2.014225e-07, 2.499519e-07, 3.485983e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 200
# )
# 
# # Convert x to factor if needed for clean spacing
# df_200_8$x <- factor(df_200_8$x)
# 
# df_200_7 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9931703, 0.9931734, 0.9931732, 0.9933937, 0.9931760,
#         0.9933236, 0.9931255, 0.9935360, 0.9931597, 0.9933033,
#         0.9932149, 0.9931922, 0.9931647, 0.9930886, 0.9931526, 
#         0.9932552, 0.9931892, 0.9929760, 0.9928216, 0.9927226),
#   sd = sqrt(c(2.387314e-07, 2.459100e-07, 1.074160e-07, 1.660072e-07, 1.627875e-07,
#               1.161549e-07, 1.199895e-07, 1.999477e-07, 1.908565e-07, 7.717771e-08,
#               2.890389e-07, 3.437699e-07, 1.346382e-07, 4.799738e-07, 2.724402e-07, 
#               1.417829e-07, 5.139579e-07, 1.875239e-07, 2.902604e-07, 3.815000e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 200
# )
# 
# # Convert x to factor if needed for clean spacing
# df_200_7$x <- factor(df_200_7$x)
# 
# df_200_5 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9930928, 0.9931029, 0.9930838, 0.9933000, 0.9931314,
#         0.9932164, 0.9930512, 0.9934245, 0.9926277, 0.9931663,
#         0.9931614, 0.9930923, 0.9930496, 0.9930262, 0.9927616, 
#         0.9931715, 0.9930942, 0.9929389, 0.9927315, 0.9925912),
#   sd = sqrt(c(2.380990e-07, 2.811883e-07, 7.804687e-08, 1.380887e-07, 1.596658e-07,
#               1.385602e-07, 1.506668e-07, 1.521153e-07, 2.008525e-07, 7.783113e-08,
#               3.410988e-07, 3.738003e-07, 1.255841e-07, 5.983880e-07, 2.891806e-07, 
#               1.426330e-07, 5.088231e-07, 1.872784e-07, 3.010474e-07, 3.985718e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 200
# )
# 
# # Convert x to factor if needed for clean spacing
# df_200_5$x <- factor(df_200_5$x)
# 
# df_200_4 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9926974, 0.9927230, 0.9926865, 0.9929020, 0.9927783,
#         0.9926913, 0.9927150, 0.9929196, 0.9926277, 0.9925258,
#         0.9926713, 0.9926604, 0.9925677, 0.9925744, 0.9921314, 
#         0.9927524, 0.9926841, 0.9924976, 0.9922833, 0.9920420),
#   sd = sqrt(c(3.313057e-07, 1.903486e-07, 6.026302e-08, 1.533625e-07, 2.072384e-07,
#               2.549396e-07, 1.360316e-07, 2.637165e-07, 2.008525e-07, 1.563458e-07,
#               3.273993e-07, 3.468608e-07, 1.115782e-07, 4.691303e-07, 5.184709e-07, 
#               2.143356e-07, 6.017132e-07, 1.459019e-07,3.930007e-07, 3.493351e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 200
# )
# 
# # Convert x to factor if needed for clean spacing
# df_200_4$x <- factor(df_200_4$x)
# 
# df_400_8 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9926275, 0.9924494, 0.9927289, 0.9924561, 0.9926964,
#         0.9927170, 0.9924875, 0.9928824, 0.9928032, 0.9924508,
#         0.9926056, 0.9925090, 0.9927561, 0.9924486, 0.9918460,
#         0.9925349, 0.9926726, 0.9922704, 0.9921367, 0.9919606),
#   sd = sqrt(c(1.749395e-07, 2.056222e-07, 7.895957e-08, 1.652997e-07, 1.361566e-07,
#               1.822957e-07, 1.482524e-07, 2.104794e-07, 1.928414e-07, 1.890905e-07,
#               1.972698e-07, 2.781873e-07, 2.328016e-07, 2.045380e-07, 1.662257e-07,
#               1.092738e-07, 1.727383e-07, 1.559446e-07, 1.731006e-07, 3.338247e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 400
# )
# 
# # Convert x to factor if needed for clean spacing
# df_400_8$x <- factor(df_400_8$x)
# 
# df_400_7 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9926239, 0.9924452, 0.9927426, 0.9924571, 0.9926914,
#         0.9927164, 0.9924690, 0.9929000, 0.9928159, 0.9923940,
#         0.9926056, 0.9925076, 0.9927402, 0.9923805, 0.9917389,
#         0.9925349, 0.9926709, 0.9922704, 0.9921367, 0.9919629),
#   sd = sqrt(c(1.750542e-07, 2.037097e-07, 8.051375e-08, 1.643939e-07, 1.390547e-07,
#               1.806037e-07, 1.489711e-07, 2.158329e-07, 1.293941e-07, 1.792717e-07,
#               1.972698e-07, 2.654226e-07, 2.216617e-07, 2.499021e-07, 1.744260e-07,
#               1.092738e-07, 1.741955e-07, 1.559446e-07, 1.731006e-07, 3.302295e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 400
# )
# 
# # Convert x to factor if needed for clean spacing
# df_400_7$x <- factor(df_400_7$x)
# 
# df_400_5 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9925579, 0.9923647, 0.9927221, 0.9923847, 0.9926359,
#         0.9926756, 0.9924087, 0.9928239, 0.9925534, 0.9921399,
#         0.9925433, 0.9924439, 0.9926899, 0.9922552, 0.9915602,
#         0.9924671, 0.9925890, 0.9922217, 0.9920516, 0.9919131),
#   sd = sqrt(c(1.557659e-07, 2.517004e-07, 8.268794e-08, 1.967741e-07, 1.661551e-07,
#               1.691365e-07, 1.797051e-07, 1.946813e-07, 2.300729e-07, 2.281220e-07,
#               1.658461e-07, 2.128074e-07, 2.152595e-07, 2.240917e-07, 1.646969e-07,
#               1.140390e-07, 1.795030e-07, 1.240123e-07, 1.862590e-07, 3.569773e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 400
# )
# 
# # Convert x to factor if needed for clean spacing
# df_400_5$x <- factor(df_400_5$x)
# 
# df_400_4 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9920932, 0.9919453, 0.9924498, 0.9920145, 0.9923265,
#         0.9923027, 0.9920341, 0.9925425, 0.9921692, 0.9914460,
#         0.9922233, 0.9920781, 0.9923284, 0.9919225, 0.9910387,
#         0.9920643, 0.9922044, 0.9917942, 0.9917035, 0.9914922),
#   sd = sqrt(c(1.776316e-07, 3.323934e-07, 7.732921e-08, 4.359305e-07, 1.977048e-07,
#               2.178169e-07, 2.268436e-07, 1.949454e-07, 2.673348e-07, 2.118785e-07,
#               1.607366e-07, 2.716196e-07, 1.786500e-07, 2.381705e-07, 2.453483e-07,
#               1.423955e-07, 2.180441e-07, 6.684048e-08, 2.491977e-07, 4.058565e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 400
# )
# 
# # Convert x to factor if needed for clean spacing
# df_400_4$x <- factor(df_400_4$x)
# 
# df_1000_8 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c( 0.9923354, 0.9923601, 0.9924817, 0.9923196, 0.9921489,
#          0.9922108, 0.9923334, 0.9920620, 0.9918252, 0.9914047,
#          0.9920139, 0.9923297, 0.9922957, 0.9918396, 0.9913304,
#          0.9922048, 0.9922374, 0.9921161, 0.9918686, 0.9915913),
#   sd = sqrt(c(6.856122e-08, 1.455216e-07, 1.400744e-07, 2.645276e-07, 2.185404e-07,
#               2.939525e-07, 8.366131e-08, 1.896804e-07, 2.726159e-07, 2.474518e-07,
#               1.851258e-07, 3.703346e-07, 2.144926e-07, 2.417547e-07, 6.785674e-07 ,
#               3.615731e-07, 2.886637e-07, 2.146313e-07, 2.576833e-07, 3.104953e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 1000
# )
# 
# # Convert x to factor if needed for clean spacing
# df_1000_8$x <- factor(df_1000_8$x)
# 
# df_1000_7 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9923352, 0.9923601, 0.9924845, 0.9923108, 0.9921527,
#         0.9922108, 0.9923212, 0.9920581, 0.9918178, 0.9913472,
#         0.9920091, 0.9923183, 0.9922957,  0.9918396,  0.9913242,
#         0.9922048, 0.9922334, 0.9921138, 0.9918691,  0.9915913),
#   sd = sqrt(c(6.8561690e-08, 1.456252e-07, 1.446672e-07, 2.674654e-07, 2.156303e-07,
#               2.939525e-07, 8.204312e-08, 1.992801e-07, 2.710210e-07, 2.357127e-07,
#               1.827716e-07, 3.573798e-07, 2.144926e-07,  2.417547e-07,  6.770498e-07,
#               3.615731e-07, 2.882895e-07, 2.163479e-07,  2.572795e-07,  3.104953e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 1000
# )
# 
# # Convert x to factor if needed for clean spacing
# df_1000_7$x <- factor(df_1000_7$x)
# 
# df_1000_5 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9922954, 0.9922884, 0.9924225, 0.9921714, 0.9920190,
#         0.9921944, 0.9922337, 0.9920456, 0.9917330, 0.9912095,
#         0.9919860, 0.9922979,  0.9922226,  0.9917912, 0.9912728,
#         0.9921525,0.9921700, 0.9920746, 0.9918021,  0.9915026),
#   sd = sqrt(c(7.518880e-08, 1.538403e-07, 1.123810e-07, 3.042485e-07, 2.516990e-07,
#               2.932652e-07, 7.498580e-08, 1.980757e-07, 2.718864e-07, 2.314191e-07,
#               1.636649e-07,  3.570028e-07, 2.007116e-07,  2.637004e-07, 6.948986e-07,
#               3.904975e-07,  3.330338e-07, 2.142572e-07,  2.548929e-07,  3.542013e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 1000
# )
# 
# # Convert x to factor if needed for clean spacing
# df_1000_5$x <- factor(df_1000_5$x)
# 
# df_1000_4 <- data.frame(
#   x = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), times = 4),
#   y = c(0.9919259, 0.9918548, 0.9920302, 0.9917288, 0.9914755,
#         0.9918497, 0.9918950, 0.9917494, 0.9913670, 0.9907395,
#         0.9916562, 0.9919009,  0.9917934, 0.9914225, 0.9907706,
#         0.9918482, 0.9917938, 0.9917215, 0.9912995, 0.9911541),
#   sd = sqrt(c(1.494389e-07, 2.286905e-07, 8.555097e-08, 2.951696e-07, 3.913785e-07,
#               3.767020e-07, 4.241307e-08, 2.025255e-07, 2.807962e-07, 2.940117e-07, 
#               1.940751e-07, 3.732257e-07,  4.199206e-07,  3.422400e-07, 6.366096e-07,
#               4.198405e-07, 2.879157e-07, 2.352077e-07,  3.109499e-07, 4.080350e-07)),
#   method = rep(c("1% Variance Explained", "5% Variance explained", "10% Variance explained", "20% Variance explained"), each = 5),
#   Individuals = 1000
# )
# 
# # Convert x to factor if needed for clean spacing
# df_1000_4$x <- factor(df_1000_4$x)


df_all <- rbind(df_200, df_400, df_1000)

ggplot(df_all, aes(x = x, y = y, color = method, group = method)) +
  geom_point(position=position_dodge(width = 0.75, preserve = "total"),size=6) +
  geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 1, size = 1, position=position_dodge(width = 0.75, preserve = "total")) +
  coord_cartesian(ylim = c(0.98, 1))  +
  facet_wrap(~ Individuals, ncol = 1, labeller = label_both) +
  labs(
    x = "Percentage of phenotypes with effect",
    y = "Average Correlation of Simulated LD vs True LD",
    color = "Percentage variance explained"
  ) +
  scale_color_manual(values = c("#4285f4","#34a853","#fbbc05","#ea4335")) +
  ggtitle("Threshold 1e-6") + 
  theme_minimal(base_size = 13) +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(r = 15)),
    axis.title.x = element_text(size = 22, margin = margin(t = 12)),
    axis.text = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
  ) 


