%Produce final Cumulative responses Ouput
%%
[PSTH] = ComputePSTHMultiD(allclusters, "Conc",4);
%%
[smoothPSTH] = SmoothPSTH(PSTH, 200, "Conc");
%%
[normalizedPSTH] = NormalizedPSTH(smoothPSTH, "Conc");
%%
[cumulativePSTH] = CumulativePSTH(normalizedPSTH, "Conc");