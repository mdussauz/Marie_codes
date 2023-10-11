%% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% fully orthogonal to all preceeding regressors in the matrix and 0 being
% fully redundant. Having fully redundant regressors in the matrix will
% break the model, so in this example those regressors are removed. In
% practice, you should understand where the redundancy is coming from and
% change your model design to avoid it in the first place!

rejIdx = false(1,size(fullR,2));
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end
% save([fPath filesep 'regData.mat'], 'fullR', 'regIdx', 'regLabels','fullQRR','-v7.3'); %save some model variables

%% check for the presence of NaN or Inf values 
TestNaN = sum(isnan(fullR));
TestInf = sum(isinf(fullR));
if TestNaN ~= 0
    fprintf('Design matrix contains NaN')
end
if TestInf ~=0
    fprintf('Design matrix contains Inf')
end