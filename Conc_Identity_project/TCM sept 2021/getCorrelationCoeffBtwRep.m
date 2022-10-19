function [CorrelationMatrix] = getCorrelationCoeffBtwRep
tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]
odormeanFR = squeeze(nanmean(x(:,:,:,10000:14000,:),4));

for whatsmell = 1:5
    for whatconc = 1:4 
        
        all_odormeanFR (:,whatconc,whatsmell,:) = odormeanFR;
    end
end 

R = corrcoef(A,B)

end 