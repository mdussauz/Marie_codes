%plot concentration invariance test

firingRatesAverage = nanmean(smoothPSTH,5);
firstbin = -10;
lastbin = 10;

color = [0, 0.447, 0.7410; 0.8500, 0.3250, 0.0980; 0.4940, 0.1840, 0.5560; 0.6350, 0.0780, 0.1840; 0.3010, 0.7450, 0.9330];

%% Invariance test plot for whole odor period
figure(1)
for whatsmell = 1:5
    FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
    subplot(1,3,1)
    X1(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:14000), 3)); % Concentration 1
    Y1(whatsmell,:) = squeeze(mean(FRAverageConc(:,2,10000:14000), 3)); % Concentration 2
    scatter(X1(whatsmell,:),Y1(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:));hold on
    
    subplot(1,3,2)
    FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
    X2(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:14000), 3)); % Concentration 1
    Y2(whatsmell,:) = squeeze(mean(FRAverageConc(:,3,10000:14000), 3)); % Concentration 3
    scatter(X2(whatsmell,:),Y2(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:));hold on
    
    subplot(1,3,3)
    X3(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:14000), 3)); % Concentration 1
    Y3(whatsmell,:) = squeeze(mean(FRAverageConc(:,4,10000:14000), 3)); % Concentration 4
    scatter(X3(whatsmell,:),Y3(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:)); hold on
end


figure(2)
for whatsmell = 1:5
    FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
    subplot(1,3,1)
    X1(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:10300), 3)); % Concentration 1
    Y1(whatsmell,:) = squeeze(mean(FRAverageConc(:,2,10000:10300), 3)); % Concentration 2
    scatter(X1(whatsmell,:),Y1(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:));hold on
    
    subplot(1,3,2)
    FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
    X2(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:10300), 3)); % Concentration 1
    Y2(whatsmell,:) = squeeze(mean(FRAverageConc(:,3,10000:10300), 3)); % Concentration 3
    scatter(X2(whatsmell,:),Y2(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:));hold on
    
    subplot(1,3,3)
    X3(whatsmell,:) = squeeze(mean(FRAverageConc(:,1,10000:10300), 3)); % Concentration 1
    Y3(whatsmell,:) = squeeze(mean(FRAverageConc(:,4,10000:10300), 3)); % Concentration 4
    scatter(X3(whatsmell,:),Y3(whatsmell,:), 25, 'filled', 'MarkerFaceColor',color(whatsmell,:)); hold on
    
end



figure(1)
title( ['Full odor period'])
set(figure(1),'Renderer','painters')

subplot(1,3,1)
xlabel('Conc 1')
ylabel('Conc 2')
x = X1(:);
y = Y1(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1, 'LineStyle', '--') %linear regression
plot([0 50], [0 50])

subplot(1,3,2)
xlabel('Conc 1')
ylabel('Conc 3')
x = X2(:);
y = Y2(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1, 'LineStyle', '--') %linear regression
plot([0 50], [0 50])

subplot(1,3,3)
xlabel('Conc 1')
ylabel('Conc 4')
x = X3(:);
y = Y3(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1,'LineStyle', '--') %linear regression
plot([0 50], [0 50])

figure(2)
title( ['300 ms of odor period'])
set(figure(2),'Renderer','painters')
subplot(1,3,1)
xlabel('Conc 1')
ylabel('Conc 2')
x = X1(:);
y = Y1(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1, 'LineStyle', '--') %linear regression
plot([0 50], [0 50])

subplot(1,3,2)
xlabel('Conc 1')
ylabel('Conc 3')
x = X2(:);
y = Y2(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1, 'LineStyle', '--') %linear regression
plot([0 50], [0 50])

subplot(1,3,3)
xlabel('Conc 1')
ylabel('Conc 4')
x = X3(:);
y = Y3(:);
b1 = x\y; % slope or regression coefficient
yCalc1 = b1*x;
plot(x,yCalc1, 'LineStyle', '--') %linear regression
plot([0 50], [0 50])



