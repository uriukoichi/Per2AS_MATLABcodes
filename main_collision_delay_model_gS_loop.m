%% Fig. 5G,H,I

global Ks Ka ns
global fS muS gS hS
global vA fA muA

mA = 1; % id for Per2AS
mS = 2; % id for Per2 
pS = 3; % id for PER2 protein

muS = log(2.0)/2.0; % degradation rate of Per2 mRNA

f = 1; % probability of collision and detachment
fS = f;

ns = 3.0; % nonlinearlity

Ks = 0.05; % ratio of RNAP detachment rate to recruitment rate for Per2

%gS = 14.0/muS; % translatin rate defaul value
gS = 14.0*1/muS; % translatin rate
hS = 1.0; % degradation rate of PER2

tau = 6.2; % time delay

vA = 1.0; % transcription rate of Per2AS
muA = 1.0; % degradation rate of Per2AS

fA = f; % probability of collision and detachment for Per2AS
Ka = 4.0; % ratio of RNAP detachment rate to recruitment rate for Per2AS


Tau = [tau, tau, tau]; % delay matrix
y0 = [0, 0, 0]; % initial condition

Tmax = 720;

ExpL = 0:0.2:3;

gSL = 14.0*10.^ExpL'/muS;
ngSL = length(gSL);

Per2_ave_data = zeros(ngSL, 1);
Per2AS_ave_data = zeros(ngSL, 1);
PER2_ave_data = zeros(ngSL, 1);
period_data = zeros(ngSL, 1);

for i = 1:ngSL
    
    gS = gSL(i);
    sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, Tmax]);

    %% period calculation
    T = 0:0.01:sol.x(end);
    Y = deval(sol, T);
    [pks, locs] = findpeaks(Y(mS,:)', T');
    nlocs = size(locs,1);
    Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
    aTp = mean(Tp);
    period_data(i) = aTp;

    %% calculation of transcriptional activities
    XS = 1./(1+Ks.*(Y(pS,:).^ns));
    XA = (1/(1+Ka))*ones(1,length(T));

    %% calculation of average level
    Tave = Tmax - 10*24;
    [mA_ave, mS_ave, pS_ave] = ave_cal(sol, Tave);
    Per2_ave_data(i) = mS_ave;
    Per2AS_ave_data(i) = mA_ave;
    PER2_ave_data(i) = pS_ave;

end

%% PLOT gS vs average value
figure
t = tiledlayout(1,3);

nexttile
semilogx(gSL, PER2_ave_data/PER2_ave_data(1), 'bo-');
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
pbaspect([1,1,1]);

nexttile
semilogx(gSL, Per2AS_ave_data/Per2AS_ave_data(1), 'ro-');
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'YLim',[0,4]);
set(gca,'YTick',0:1:4);
% set(gca,'XLim',[10,10^4]);
pbaspect([1,1,1]);

nexttile
semilogx(gSL, Per2_ave_data/Per2_ave_data(1), 'ko-');
set(gca,'FontName','Arial');
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'YLim',[0,4]);
set(gca,'YTick',0:1:4);
pbaspect([1,1,1]);

t.TileSpacing = 'compact';
t.Padding = 'compact';


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



