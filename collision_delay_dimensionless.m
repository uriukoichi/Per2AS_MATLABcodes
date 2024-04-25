function dy = collision_delay_dimensionless(t, y, Z)

% function to define delay differential equations
% including collision
% quasi steady sate approximation

global Ks Ka ns
global fS muS gS hS
global vA fA muA

mA = 1; % id for Per2AS
mS = 2; % id for Per2 
pS = 3; % id for PER2 protein

dy = zeros(3,1);

ylag = Z(:,pS); % delay for PER2 protein translation
mSd = ylag(mS);

% transcriptional activities
Xs = 1.0/(1.0 + Ks*(y(pS)^ns));
Xa = 1.0/(1.0 + Ka);

% delay differential equation
dy(mS) = muS*( Xs*(1.0-Xa) + (1.0-fS)*Xs*Xa - y(mS) );
dy(mA) = muS*( vA*Xa*(1.0-Xs) + vA*(1.0-fA)*Xa*Xs - muA*y(mA) );
dy(pS) = muS*( gS*mSd - hS*y(pS) );

end