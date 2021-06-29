%% This is to plot the dPCs directly
% 


T = 20000;
%T = 4000; % when considering only odor timewindow
ncomp = 20;
ToSmooth = 0;
sp = 25;

dpc = W(:,1:ncomp)'*X;
x = reshape(dpc,ncomp, 5, 4, T);

%which_comps = [4 15 18]; % for APC identity start trial/odor end odor
which_comps = [3 8 17]; % for APC conc
%which_comps = [6 10 12]; % for APC identity
%which_comps = [4 5 9]; % for AON identity
%which_comps = [2 8 15]; % for AON concentration
%which_comps = [3 11 20]; % Change these numbers to the top 3 components you need to plot;
% this will depend on whether you are plotting identity or concentration
% subspaces

color = [165,0,38;253,174,97;116,173,209;69,117,180;49,54,149]/255;
lw = ([0.5 1 2 4]);
a = lw;
view(-25,10) % for turning the viewing angle of 3D plots

test = 9000:14000; %1s before odor start to odor end
%test = : ; % initial value

for odor = [1 2 3 4 5]
        hold on;
        for dil = 1:4
            temp2 = squeeze(x(which_comps,odor,dil,test));
            %temp2 = squeeze(x(which_comps,odor,dil,:));
            if (ToSmooth == 1)
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
                set(h1,'Color',color(odor,:),'LineWidth',a(dil));
            else
                if size(which_comps,2)==3
                    h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                elseif size(which_comps,2)==2
                   h1 = plot(squeeze(temp2(1,:)),squeeze(temp2(2,:))); %if only 2 PC to plot
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                end 
            end
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            if size(which_comps,2)==3
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
            end 
        end
end
