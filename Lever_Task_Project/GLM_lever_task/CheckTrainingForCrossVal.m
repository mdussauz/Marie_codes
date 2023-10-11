%CheckingTrainingForCrossVal 
ridgeFolds = 10;
Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
randIdx = randperm(size(spikeTrace,1)); %generate randum number index
foldCnt = floor(size(spikeTrace,1) / ridgeFolds);



 for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            test = ((iFolds - 1)*foldCnt) + (1:foldCnt); 
            test2 = dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt)));
            y_here = ones(size(test2));
            y_here2 = ones(size(test));
            figure(1)
            plot(test2, y_here* iFolds); hold on; 
            figure(2)
            plot(test, y_here2* iFolds); hold on; 
 end 
 
%  test = ((9 - 1)*foldCnt) + (1:foldCnt); 
%  y_here = ones(size(test));
%  plot(test7, y_here); hold on; 
%  plot(test8, y_here*2); hold on; 
%  plot(test9, y_here*3)

figure(3)
plot(randIdx, ones(size(randIdx)))