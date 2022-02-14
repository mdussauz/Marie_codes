%plot FR average for each unit
units_per_fig = 5;
x = -10:1.0000e-03:10-1.0000e-03;
Nneurons = size(smoothPSTH, 1);
firingRatesAverage = nanmean(smoothPSTH,5);

% create a default color map ranging from red to light pink
length = 5;
red = [0, 0, 1];
color_chosen = [0, 1, 1];
colors_p = [linspace(red(1),color_chosen(1),length)', linspace(red(2),color_chosen(2),length)', linspace(red(3),color_chosen(3),length)'];

i = 1; j =1; All_max_AV =[]; All_min_AV =[];% initialize

for n = 1:Nneurons % for every cell 
    for whatsmell = 1:5
        figure(j);
        subplot(units_per_fig,5,i);
        for whatconc = 1:4
            FRAverage = squeeze(firingRatesAverage(n,whatsmell,whatconc,:));
            p = plot(x,FRAverage, 'Color',colors_p(5-whatconc,:)); hold on 
            
            max_AV = max(FRAverage);
            All_max_AV = [All_max_AV, max_AV];
            min_AV = min(FRAverage);
            All_min_AV = [All_min_AV, min_AV];
            
            xlim([-9.8 9.8])
            ylim([min(All_min_AV)-0.5 max(All_max_AV)+0.5])
            p.LineWidth = 1;
        end
                i = i +1;
                All_max_AV =[]; 
                All_min_AV =[];
    end
                if mod(n, units_per_fig) == 0
                i=1;
                j = j +1;
            end  
end
