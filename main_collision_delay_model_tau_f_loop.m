%% Code for Fig. 2FG 

tic

global Ks Ka ns
global fS muS gS hS
global vA fA muA

mA = 1;
mS = 2;
pS = 3;

muS = log(2.0)/2.0;

ns = 3.0;

Ks = 0.05;

gS = 14.0/muS;
hS = 1.0;

fS = 0.0;
fA = 0.0;

tau = 6.2;

vA = 1.0;
muA = 1.0;

Ka = 4.0;

q_amp = 0.01; % threshold for amplitude

y0 = [0, 0, 0];

taul = 2:0.1:10; % tau list
fl = 0:0.1:1.0; % f list

AL_Per2 = zeros(length(taul),length(fl)); % amplitude list for Per2
AL_Per2AS = zeros(length(taul),length(fl)); % amplitude list for Per2AS
AL_PER2 = zeros(length(taul),length(fl)); % amplitude list for PER2 protein
PL = zeros(length(taul),length(fl)); % period list

for i = 1:length(taul)
    for j = 1:length(fl)

    tau = taul(i);
    Tau = [tau, tau, tau];

    fS = fl(j);
    fA = fl(j);

    sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, 2400]);

    %% amplitude & period calculation
    T = 0:0.01:sol.x(end);
    Y = deval(sol, T);
    
    %% Per2AS peak and trough
    [pks,~] = findpeaks(Y(mA,:)', T');
    [trs,~] = findpeaks(-Y(mA,:)', T');
    npks = size(pks,1);
    apk = mean(pks(end-round(npks/3)+1:end));
    ntrs = size(trs,1);
    atr = -mean(trs(end-round(ntrs/3)+1:end));
    
    AL_Per2AS(i,j) = apk - atr;

    % no oscillation of Per2AS for f = 0
    if j == 1
        AL_Per2AS(i,j) = 0.0;
    end

    %% Per2 peak and trough
    [pks, locs] = findpeaks(Y(mS,:)', T');
    [trs,~] = findpeaks(-Y(mS,:)', T');
    npks = size(pks,1);
    apk = mean(pks(end-round(npks/3)+1:end));
    ntrs = size(trs,1);
    atr = -mean(trs(end-round(ntrs/3)+1:end));

    AL_Per2(i,j) = apk - atr;
    
    if AL_Per2(i,j) > q_amp
        % Period calculation
        nlocs = size(locs,1);
        Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
        aTp = mean(Tp);

        PL(i,j) = aTp;
    else
        PL(i,j) = NaN; % if there is no oscillation, we set NaN
    end

    %% PER2 protein peak and trough
    [pks, ~] = findpeaks(Y(pS,:)', T');
    [trs,~] = findpeaks(-Y(pS,:)', T');
    npks = size(pks,1);
    apk = mean(pks(end-round(npks/3)+1:end));
    ntrs = size(trs,1);
    atr = -mean(trs(end-round(ntrs/3)+1:end));

    AL_PER2(i,j) = apk - atr;
    

    end
end

%% Plotting
t = tiledlayout(1,4);
ax1 = nexttile;
[X, Y] = meshgrid(fl, taul);
contourf(X, Y, AL_Per2);
colormap(ax1,"default");
clim([0, 1])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[2, 10]);
set(gca,'YTick',2:1:10);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10, 'Ticks', [0.0,0.5,1],'TickDirection','out', 'Color','k')

hold on
yline(6.2);
hold off

ax2 = nexttile;
contourf(X, Y, AL_Per2AS);
colormap(ax2,"parula");
clim([0, 0.2])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[2, 10]);
set(gca,'YTick',2:1:10);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'Ticks', [0.0,0.1,0.2],'TickDirection','out','Color','k')

hold on
yline(6.2);
hold off

ax3 = nexttile;
contourf(X, Y, AL_PER2);
colormap(ax3,"parula");
clim([0, 32.0])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[2, 10]);
set(gca,'YTick',2:1:10);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'Ticks', [0.0,16,32],'TickDirection','out','Color','k')
%colorbar('northoutside','FontName', 'Arial', 'FontSize',10, 'TickDirection','out', 'Color','k')

hold on
yline(6.2);
hold off


ax4 = nexttile;
contourf(X, Y, PL);
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
colormap(ax4, mycolormap);
clim([16, 32])
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'XTick',0:0.2:1);
set(gca,'YLim',[2, 10]);
set(gca,'YTick',2:1:10);
pbaspect([1,1,1]);
colorbar('northoutside','FontName', 'Arial', 'FontSize',10,'Ticks', [16,20,24,28,32],'TickDirection','out','Color','k')

hold on
yline(6.2);
hold off

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% OUTPUT
% Str = strcat('tau_vs_mA_mS_peak_trough_f_',num2str(fS,'%.2f'),'.dat');
% fileID = fopen(Str, 'w');
% fprintf(fileID, '%6.4f %6.4f %6.4f %6.4f %6.4f\n', [taul; PT']);
% fclose(fileID);
% 
% Str = strcat('tau_vs_Per2_period.dat');
% fileID = fopen(Str, 'w');
% fprintf(fileID, '%6.4f %6.4f\n', [taul; PL']);
% fclose(fileID);
% 

save("f_vs_tau_vs_amplitude_period.mat")

toc


