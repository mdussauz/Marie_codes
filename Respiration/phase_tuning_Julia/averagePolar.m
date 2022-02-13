function[averageTheta, averageRho] = averagePolar(polars)
% Written by JW
% Get the vector average of a list of polar phase vectors
% Input:
%   -polars = list of thetas
    [x, y] = pol2cart(polars, ones('like', polars));
    sumX = sum(x);
    sumY = sum(y);
    [averageTheta, rho] = cart2pol(sumX, sumY);
    averageRho = rho/length(polars);
end