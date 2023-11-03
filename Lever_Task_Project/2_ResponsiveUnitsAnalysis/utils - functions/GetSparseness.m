function [LifetimeSparsness, PopulationSparsness] = GetSparseness(SparseVar)
Num1 = (nansum(bsxfun(@rdivide,SparseVar, sum(~isnan(SparseVar))))).^2;
Num2 = nansum(bsxfun(@rdivide,SparseVar.^2,sum(~isnan(SparseVar))));
SparseNumerator = 1-Num1./Num2;
SparseDenominator = 1-(1./sum(~isnan(SparseVar)));
LifetimeSparsness = squeeze(SparseNumerator./SparseDenominator);

SparseVar = permute(SparseVar,[2,1,3]);
Num1 = (nansum(bsxfun(@rdivide,SparseVar, sum(~isnan(SparseVar))))).^2;
Num2 = nansum(bsxfun(@rdivide,SparseVar.^2,sum(~isnan(SparseVar))));
SparseNumerator = 1-Num1./Num2;
SparseDenominator = 1-(1./sum(~isnan(SparseVar)));
PopulationSparsness = squeeze(SparseNumerator./SparseDenominator);
end