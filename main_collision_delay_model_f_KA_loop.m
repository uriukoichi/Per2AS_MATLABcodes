global Ks Ka ns
global fS muS gS hS
global vA fA muA

mA = 1; % id for Per2AS
mS = 2; % id for Per2 
pS = 3; % id for PER2 protein

%% Parameters for Per2 and PER2
muS = log(2.0)/2.0; % degradation rate of Per2 mRNA

ns = 3.0; % nonlinearlity

Ks = 0.05; % ratio of RNAP detachment rate to recruitment rate for Per2

gS = 14.0/muS; % translation rate of PER2 protein
hS = 1.0; % degradation rate of PER2 protein

tau = 6.2; % time delay in PER2 protein production

%% Parameters for Per2AS
vA = 1.0; % transcription rate of Per2AS
muA = 1.0; % degradation rate of Per2AS

%% Setting delay and initial condition
Tau = [tau, tau, tau];
y0 = [0, 0, 0];

Tmax = 720; % claculation time

f_list = 0:0.1:1.0; % list of probability of collision and detachment
nfL = size(f_list, 2);

Ka_list = 0.6:0.34:4.0; % list of ratio of RNAP detachment rate to recruitment rate for Per2AS
nKaL = size(Ka_list, 2);

Per2_ave_data = zeros(nfL, nKaL); % amplitude matrix of Per2
Per2AS_ave_data = zeros(nfL, nKaL); % amplitude matrix of Per2AS
PER2_ave_data = zeros(nfL, nKaL); % amplitude matrix of PER2 protein
period_data = zeros(nfL, nKaL); % period matrix

for i = 1:nfL
    for j = 1:nKaL

        fS = f_list(i); % probability of collision and detachment
        fA = f_list(i);
        Ka = Ka_list(j); % ratio of RNAP detachment rate to recruitment rate for Per2AS

        sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, Tmax]);

        %% calculation of period
        aTp = period_calculation(sol);
        period_data(j, i) = aTp;
        
        %% calculation of average level
        Tave = Tmax - 10*aTp;
        [mA_ave, mS_ave, pS_ave] = ave_cal(sol, Tave);
        Per2_ave_data(j,i) = mS_ave;
        Per2AS_ave_data(j,i) = mA_ave;
        PER2_ave_data(j,i) = pS_ave;

    end
end

%% Plotting
t = tiledlayout(1,3);
nexttile
[X, Y] = meshgrid(f_list, Ka_list);
contourf(X, Y, Per2_ave_data);
clim([0.1, 0.24])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[0.6, 4]);
set(gca,'YTick',0.6:0.68:4);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10, 'Ticks', [0.1,0.135,0.17,0.205,0.24],'TickDirection','out', 'Color','k')

nexttile
contourf(X, Y, Per2AS_ave_data);
clim([0, 0.6])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[0.6, 4]);
set(gca,'YTick',0.6:0.68:4);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'TickDirection','out','Color','k')

% nexttile
% contourf(X, Y, PER2_ave_data);
% clim([5, 9.5])
% set(gca,'FontName','Arial');
% set(gca,'XColor','k');
% set(gca,'YColor','k');
% set(gca,'XTick',0:0.2:1);
% set(gca,'YLim',[0.6, 4]);
% set(gca,'YTick',0.6:0.68:4);
% pbaspect([1,1,1]);
% colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'TickDirection','out','Color','k')

nexttile
contourf(X, Y, period_data);
clim([22, 24])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[0.6, 4]);
set(gca,'YTick',0.6:0.68:4);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'TickDirection','out','Color','k')

t.TileSpacing = 'compact';
t.Padding = 'compact';


%% OutPut
%save('fS_fA_vs_amp_period.mat');

%% Functions
% Period calculation
function aTp = period_calculation(sol)

mS = 2;

T = 0:0.1:sol.x(end);
Y = deval(sol, T);

[~, locs] = findpeaks(Y(mS,:)', T');

nlocs = size(locs,1);
Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
aTp = mean(Tp);

end

%% Average calculation
function [mA_ave, mS_ave, pS_ave] = ave_cal(sol, Tave)

mA = 1;
mS = 2;
pS = 3;

T = 0:0.01:sol.x(end);
Y = deval(sol, T);

I = find(T(:) >= Tave);

mA_ave = mean(Y(mA, I));
mS_ave = mean(Y(mS, I));
pS_ave = mean(Y(pS, I));

end



