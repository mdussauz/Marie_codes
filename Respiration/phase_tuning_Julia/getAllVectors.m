function[airthetas, airrhos, odorthetas, odorrhos] = getAllVectors(cluster, phaseStarts)
% Gets all the average air and odor vectors (in polar form) for a cluster
% Input:
% -cluster
% -phaseStarts: a list of respiration phase starts
airthetas = zeros(20, 1, 'double');
airrhos = zeros(20, 1, 'double');
odorthetas = zeros(20, 1, 'double');
odorrhos = zeros(20, 1, 'double');
for ind = 1:20
    [airpolars, odorpolars] = getAirandOdorPolars(cluster, ind, phaseStarts);
    [airTheta, airRho] = averagePolar(airpolars);
    [odorTheta, odorRho] = averagePolar(odorpolars);
    airthetas(ind) = airTheta;
    airrhos(ind) = airRho;
    odorthetas(ind) = odorTheta;
    odorrhos(ind) = odorRho;
end
airthetas = reshape(airthetas, 5, 4);
airrhos = reshape(airrhos, 5, 4);
odorthetas = reshape(odorthetas, 5, 4);
odorrhos = reshape(odorrhos, 5, 4);

end