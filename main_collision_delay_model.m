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

gS = 14.0/muS; % translatin rate 
hS = 1.0; % degradation rate of PER2

tau = 6.2; % time delay

vA = 1.0; % transcription rate of Per2AS
muA = 1.0; % degradation rate of Per2AS

fA = f; % probability of collision and detachment for Per2AS
Ka = 4.0; % ratio of RNAP detachment rate to recruitment rate for Per2AS

Tau = [tau, tau, tau]; % delay matrix
y0 = [0, 0, 0]; % initial condition 

sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, 720]);

%% period calculation
T = 0:0.01:sol.x(end);
Y = deval(sol, T);
[pks, locs] = findpeaks(Y(mS,:)', T');
nlocs = size(locs,1);
Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
aTp = mean(Tp);

%% calculation of transcriptional activities
XS = 1./(1+Ks.*(Y(pS,:).^ns));
XA = (1/(1+Ka))*ones(1,length(T));

%% PLOT

Xmax = 120;

figure
t = tiledlayout(3,1);

nexttile
plot(T, XA,'r-', T, XS,'k-');
set(gca,'YLim',[0,1.05]);
set(gca,'YTick',0:0.2:1);
set(gca,'XLim',[0,Xmax]);
set(gca,'XTick',0:24:Xmax);
set(gca,'FontName','Arial')

nexttile
plot(sol.x, sol.y(mA,:),'r-', sol.x, sol.y(mS,:),'k-');
set(gca,'YLim',[0,1]);
set(gca,'XLim',[0,Xmax]);
set(gca,'XTick',0:24:Xmax);
set(gca,'YTick',0:0.2:1);
set(gca,'FontName','Arial')

nexttile
plot(sol.x, sol.y(pS,:),'b-');
set(gca,'YLim',[0,30]);
set(gca,'XLim',[0,Xmax]);
set(gca,'XTick',0:24:Xmax);
set(gca,'YTick',0:5:30);
set(gca,'FontName','Arial')

t.TileSpacing = 'compact';
t.Padding = 'compact';

%% OUTPUT
% Str = strcat('t_vs_xA_xS_mA_mS_pS_f_',num2str(f,'%.2f'),'.dat');
% fileID = fopen(Str, 'w');
% fprintf(fileID, '%4.2f %6.4f %6.4f %6.4f %6.4f %6.4f\n', [T; XA; XS; Y]);
% fclose(fileID);




