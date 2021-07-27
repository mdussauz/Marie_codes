%Produce final Cumulative responses Ouput
%%
[PSTH] = ComputePSTHMultiD(allclusters, "Conc",4);
%%
[smoothPSTH] = SmoothPSTH(PSTH, 200, "Conc");
%%
[cumulativePSTH] = CumulativePSTH(smoothPSTH, "Conc");