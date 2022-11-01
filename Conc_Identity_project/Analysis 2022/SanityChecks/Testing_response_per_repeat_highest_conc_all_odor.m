%Test plot PSTH response

conc = 4;
%% raster
for cluster = 1: 36
    count = 1;
    for odor = 1:5
        for rep = 1:7
            figure(cluster)
            subplot(5,7,count)
            plot(squeeze(PSTH(cluster, odor, conc, :, rep)))
            count = count +1
            
        end
    end
end

%% smooth PSTH
for cluster = 1: 36
    count = 1;
    for odor = 1:5
        for rep = 1:7
            figure(cluster)
            subplot(5,7,count)
            plot(squeeze(smoothPSTH(cluster, odor, conc, :, rep)))
            count = count +1
            
        end
    end
end
