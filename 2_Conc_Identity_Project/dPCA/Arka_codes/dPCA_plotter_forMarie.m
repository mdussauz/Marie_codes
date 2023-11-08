%% Plotting the dPCs - each odor is a color and concentartions are line thickness

T = 500; % Change this to reflect the length of time axis
ncomp = 20;
ToSmooth = 1;
odor_p = 1:400; % Change this to indicate the odor period
sp = 25;
X = firingRatesAverage(:,:);


dpc = W(:,1:ncomp)'*Xtest;
x = reshape(dpc,ncomp, 5, 4, T);
which_comps = [3 11 20]; % Change these numbers to the top 3 components you need to plot;
% this will depend on whether you are plotting identity or concentration
% subspaces


figure(103);subplot(122);hold on;axis('square');
% color = 'bgcrk';
color = [1 0 0; 0 1 0; 0 0 1; 0.9 0.9 0.9; 1 0 1];
a = [1 1.2 1.7 2.5];
lw = ([0.5 1 2 4]);
view(-25,10) % for turning the viewing angle of 3D plots

for odor = [1 2 3 4 5]
        hold on;
        for dil = 1:4
            temp2 = squeeze(x(which_comps,odor,dil,:));
            if (ToSmooth == 1)
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
%                 set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
                set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
            else
                h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
%                 set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
                set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
            end
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
        end
end


%% This is to plot the dPCs directly

T = 500;
figure;
DPC = reshape(dpc,[20 5 4 T]);
color = 'bgcrk';
which_comps = [3 11 20]; % Change these numbers to the top 3 components you need to plot;
% this will depend on whether you are plotting identity or concentration
% subspaces

subplot(2,1,2)
for odor = 1:5
    for dil = 1:4
        temp = smooth(squeeze(DPC(which_comps,odor,dil,:)),10);
        plot(temp,'Color',color(odor),'LineWidth',dil/2);
        hold on;
    end
end