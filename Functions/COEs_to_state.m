%% Convert COEs to state vector
% Takes in 6 COEs and returns corresponding r and v vector

function [r, v] = COEs_to_state(COEs, mu)

a = COEs(1);
e = COEs(2);
inc = COEs(3);
RAAN = COEs(4);
omega = COEs(5);
theta = COEs(6);

h = sqrt(a*mu*(1-e^2));

rp = (h^2/mu) * (1/(1 + e*cosd(theta))) * (cosd(theta)*[1;0;0] + sind(theta)*[0;1;0]);
vp = (mu/h) * (-sind(theta)*[1;0;0] + (e + cosd(theta))*[0;1;0]);


Q = Qfunc(omega, inc, RAAN);
Q = Q';

r = (Q*rp)';
v = (Q*vp)';

end

function Q = Qfunc(omega, inc, RAAN)

R3_omega = [cosd(omega) sind(omega) 0; -sind(omega) cosd(omega) 0; 0 0 1];
R1_inc = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
R3_RAAN = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];

Q = R3_omega * R1_inc * R3_RAAN;

end