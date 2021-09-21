function [CorrelationMatrix] = getCorrelationCoeffBtwRep
tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]

for whatsmell = 1:5
    for whatconc = 1:4 
        for rep = 1:5
        meanOdorFR = squeeze(mean(x(:,whatsmell,whatconc,10000:14000,rep),4  ));


R = corrcoef(A,B)

end 