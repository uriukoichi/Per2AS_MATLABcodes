global Ks Ka ns
global fS muS gS hS
global vA fA muA

mA = 1;
mS = 2;
pS = 3;

%% Parameters for Per2 and PER2
muS = log(2.0)/2.0; % degradation rate of Per2 mRNA
%muS = log(2.0)/0.4; % for Knock Down of Per2

f = 1.0;

fS = f; % strength of interference on Per2 by Per2AS

ns = 3.0; % Hill exponent

Ks = 0.05; %

gS = 14.0/muS; % translation rate of PER2 protein
hS = 1.0; % degradation rate of PER2 protein

tau = 6.2; % time delay in PER2 protein production

%% Parameters for Per2AS
vA = 1.0; % transcription rate of Per2AS
fA = f; % strength of interference on Per2AS by Per2
muA = 1.0; % degradation rate of Per2AS
Ka = 4.0; % ratio between kon and koff for Per2AS

%% Setting delay and initial condition
Tau = [tau, tau, tau];
y0 = [0, 0, 0];

Tmax = 2400;

f_list = 0:0.05:1;
f_list = f_list';
nfL = size(f_list, 1);

Ka_list = 0.0:0.03:6.0;
Ka_list = Ka_list';
nKaL = size(Ka_list, 1);

Per2_ave_data = zeros(nKaL,nfL); % average Per2
Per2_amp_data = zeros(nKaL,nfL); % Per2 amplitude
Per2AS_ave_data = zeros(nKaL,nfL); % average Per2AS
PER2_ave_data = zeros(nKaL,nfL); % average PER2 protein
period_data = zeros(nKaL,nfL); % period

%%

for i = 1:nKaL
    for j = 1:nfL

    Ka = Ka_list(i);
    fS = f_list(j);
    fA = f_list(j);

    sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, Tmax]);

    %% calculation of period
    [aTp, apk_Per2, atr_Per2] = period_calculation(sol);
    Per2_amp_data(i,j) = apk_Per2-atr_Per2;

    if Per2_amp_data(i,j) < 0.001
        aTp = NaN;
    end

    period_data(i,j) = aTp;

    %% calculation of average level
    if ~isnan(aTp)
        Tave = Tmax - 10*aTp;
    else
        Tave = Tmax - 10*24;
    end

    [mA_ave, mS_ave, pS_ave] = ave_cal(sol, Tave);
    Per2_ave_data(i,j) = mS_ave;
    Per2AS_ave_data(i,j) = mA_ave;
    PER2_ave_data(i,j) = pS_ave;

    end
end

%% Plotting

myyrange = [0.0, 6];
myytick = 0.0:1:6;

t = tiledlayout(2,2);
ax1 = nexttile;
[X, Y] = meshgrid(f_list, Ka_list);
contourf(X, Y, Per2_amp_data);
clim([0.0, 0.7])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',myyrange);
set(gca,'YTick',myytick);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10, 'Ticks', [0.0,0.35,0.7],'TickDirection','out', 'Color','k')

hold on
yline(4);
yline(0.67);
hold off

ax2 = nexttile;
contourf(X, Y, Per2_ave_data);
clim([0.1, 0.24])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',myyrange);
set(gca,'YTick',myytick);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10, 'Ticks', [0.1,0.17,0.24],'TickDirection','out', 'Color','k')

hold on
yline(4);
yline(0.67);
hold off

ax3 = nexttile;
contourf(X, Y, Per2AS_ave_data);
clim([0, 1.0])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',myyrange);
set(gca,'YTick',myytick);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'Ticks', [0.0,0.5,1.0],'TickDirection','out','Color','k')

hold on
yline(4);
yline(0.67);
hold off

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

ax4 = nexttile;
contourf(X, Y, period_data);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
colormap(ax4, mycolormap);
clim([22, 26])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',myyrange);
set(gca,'YTick',myytick);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'TickDirection','out','Color','k')

hold on
yline(4);
yline(0.67);
hold off

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% OutPut
%save('f_Ka_vs_amp_ave_period.mat');

%% Functions
% Period calculation
function [aTp, apk_Per2, atr_Per2] = period_calculation(sol)

mS = 2;

T = 0:0.01:sol.x(end);
Y = deval(sol, T);

[pks, locs] = findpeaks(Y(mS,:)', T');
[trs, ~] = findpeaks(-Y(mS,:)', T');

npks = size(pks,1);
if npks > 0
    apk_Per2 = mean(pks(end-round(npks/3)+1:end));
    ntrs = size(trs,1);
    atr_Per2 = -mean(trs(end-round(ntrs/3)+1:end));

    nlocs = size(locs,1);
    Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
    aTp = mean(Tp);
else
    apk_Per2 = Y(mS, end);
    atr_Per2 = Y(mS, end);
    aTp = NaN;
end

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



