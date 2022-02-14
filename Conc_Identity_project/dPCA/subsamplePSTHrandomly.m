function [subsampledPSTH] = subsamplePSTHrandomly(PSTH,neurons)

%subsample neurons to be included in PSTH
%neurons = number of neurons to include in new PSTH

Nneurons= size(PSTH,1);

whichneurons = datasample([1:Nneurons],neurons,'Replace',false);
if ndims(PSTH) == 4
    subsampledPSTH = PSTH(whichneurons, :,:,:);
elseif ndims(PSTH) == 5
    subsampledPSTH = PSTH(whichneurons, :,:,:,:);
end

end