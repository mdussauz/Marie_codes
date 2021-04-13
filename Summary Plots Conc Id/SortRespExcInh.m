% Sort neurons into responsive excitatory and inhibitory neurons

% ttest (x, y, , 'alpha', 0.05,'tail', 'right') % right for mean of x is greater than mean of y? 
firingRatesAverage = nanmean(smoothPSTH,5);
Nneurons = size(smoothPSTH, 1);


for n = 1:Nneurons %for each neuron
% airmean = squeeze(nanmean(firingRatesAverage(n,:,:,6000:10000),4));
% odormean = squeeze(nanmean(firingRatesAverage(n,:,:,10000:14000),4));
% t_test_exc(n,:,:) = ttest2 (odormean, airmean, 'alpha', 0.05,  'Tail', 'right');
% t_test_inh(n,:,:) = ttest2 (odormean, airmean, 'alpha', 0.05,  'Tail', 'left');
air = squeeze(firingRatesAverage(n,:,:,6000:10000));
odor = squeeze(firingRatesAverage(n,:,:,10000:14000));
t_test_exc(n,:,:) = ttest2 (odor, air, 'alpha', 0.05,'Dim', 3,  'Tail', 'right');
t_test_inh(n,:,:) = ttest2 (odor, air, 'alpha', 0.05,'Dim',3,  'Tail', 'left');
end

IdxExcNeuron = find(t_test_exc(:,:,4)==1);
IdxInhNeuron = find(t_test_inh(:,:,4)==1);

for n = 1:Nneurons
    for whatsmell = 1:5
        for whatconc = 1:4
            
        end
    end
end