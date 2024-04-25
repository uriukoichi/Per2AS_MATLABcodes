%% Code for Fig. 2FG 

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

tau = 6.2;

vA = 1.0;
muA = 1.0;

Ka = 4.0;

Tau = [tau, tau, tau];
y0 = [0, 0, 0];

fl = 0:0.1:1;
PT = zeros(length(fl),4); % peak and trough
PL = zeros(length(fl),1); % period list

for i = 1:length(fl)

    fS = fl(i);
    fA = fl(i);

    sol = dde23(@collision_delay_dimensionless, Tau, y0, [0, 720]);

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
    
    PT(i,1) = apk;
    PT(i,2) = atr;

    % no oscillation of Per2AS for f = 0
    if i == 1
        PT(i,1) = Y(mA,end);
        PT(i,2) = Y(mA,end);
    end

    %% Per2 peak and trough
    [pks, locs] = findpeaks(Y(mS,:)', T');
    [trs,~] = findpeaks(-Y(mS,:)', T');
    npks = size(pks,1);
    apk = mean(pks(end-round(npks/3)+1:end));
    ntrs = size(trs,1);
    atr = -mean(trs(end-round(ntrs/3)+1:end));

    % Period calculation
    nlocs = size(locs,1);
    Tp = locs(end-round(nlocs/3)+1:end)-locs(end-round(nlocs/3):end-1);
    aTp = mean(Tp);

    PL(i) = aTp;
    
    PT(i,3) = apk;
    PT(i,4) = atr;

end

%%

figure
plot(fl, PT(:,1),'r-',fl, PT(:,2),'k-');
set(gca,'FontName','Arial')

figure
plot(fl, PT(:,3),'b-',fl, PT(:,4),'k-');
set(gca,'FontName','Arial')

figure
plot(fl, PL,'b-');
set(gca,'FontName','Arial')

%% OUTPUT
% Str = strcat('f_vs_mA_mS_peak_trough.dat');
% fileID = fopen(Str, 'w');
% fprintf(fileID, '%6.4f %6.4f %6.4f %6.4f %6.4f\n', [fl; PT']);
% fclose(fileID);

% Str = strcat('f_vs_Per2_period.dat');
% fileID = fopen(Str, 'w');
% fprintf(fileID, '%6.4f %6.4f\n', [fl; PL']);
% fclose(fileID);



