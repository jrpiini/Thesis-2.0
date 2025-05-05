function dstate = cowell(time, state)

% YN: 1 means yes, model that perturbation. 0 means no, ignore it

mu = 398600.44189;
J2 = 1.08262668355e-3;
J3 = -2.53265648533e-6;
J4 = -1.61962159137e-6;
J5 = -2.27296082869e-7;
J6 = 5.40681239107e-7;

r_earth = 6378.137;


x = state(1);
y = state(2);
z = state(3);

dx = state(4);
dy = state(5);
dz = state(6);

r = norm([x y z]);

ddx = -mu*x/r^3;
ddy = -mu*y/r^3;
ddz = -mu*z/r^3;


% J2
a_J2_I = (-3*J2*mu*r_earth^2*x)/(2*r^5) * (1 - ((5*z^2)/r^2));
a_J2_J = (-3*J2*mu*r_earth^2*y)/(2*r^5) * (1 - ((5*z^2)/r^2));
a_J2_K = (-3*J2*mu*r_earth^2*z)/(2*r^5) * (3 - ((5*z^2)/r^2));
a_J2 = [a_J2_I a_J2_J a_J2_K];


% J3
a_J3_I = (-5*J3*mu*r_earth^3*x)/(2*r^7) * ((3*z) - ((7*z^3)/r^2));
a_J3_J = (-5*J3*mu*r_earth^3*y)/(2*r^7) * ((3*z) - ((7*z^3)/r^2));
a_J3_K = (-5*J3*mu*r_earth^3)/(2*r^7) * ((6*z^2) - ((7*z^4)/r^2) - (3*r^2/5));
a_J3 = [a_J3_I a_J3_J a_J3_K];


% J4
a_J4_I = ((15*J4*mu*r_earth^4*x)/(8*r^7)) * (1 - (14*z^2)/(r^2) + (21*z^4)/(r^4));
a_J4_J = ((15*J4*mu*r_earth^4*y)/(8*r^7)) * (1 - (14*z^2)/(r^2) + (21*z^4)/(r^4));
a_J4_K = ((15*J4*mu*r_earth^4*z)/(8*r^7)) * (5 - (70*z^2)/(3*r^2) + (21*z^4)/(r^4));
a_J4 = [a_J4_I a_J4_J a_J4_K];

% J5
a_J5_I = ((3*J5*mu*r_earth^5*x*z)/(8*r^9)) * (35 - 210*((z^2)/(r^2)) + 231*((z^4)/(r^4)));
a_J5_J = ((3*J5*mu*r_earth^5*y*z)/(8*r^9)) * (35 - 210*((z^2)/(r^2)) + 231*((z^4)/(r^4)));
a_J5_K = (((3*J5*mu*r_earth^5*z^2)/(8*r^9)) * (105 - 315*((z^2)/(r^2)) + 231*((z^4)/(r^4)))) - ((15*J5*mu*r_earth^5)/(8*r^7));
a_J5 = [a_J5_I a_J5_J a_J5_K];

% J6
a_J6_I = ((-J6*mu*r_earth^6*x)/(16*r^9)) * (35 - 945*((z^2)/(r^2)) + 3465*((z^4)/(r^4)) - 3003*((z^6)/(r^6)));
a_J6_J = ((-J6*mu*r_earth^6*y)/(16*r^9)) * (35 - 945*((z^2)/(r^2)) + 3465*((z^4)/(r^4)) - 3003*((z^6)/(r^6)));
a_J6_K = ((-J6*mu*r_earth^6*z)/(16*r^9)) * (245 - 2205*((z^2)/(r^2)) + 4851*((z^4)/(r^4)) - 3003*((z^6)/(r^6)));
a_J6 = [a_J6_I a_J6_J a_J6_K];


a = [ddx ddy ddz] + a_J2 + a_J3 + a_J4 + a_J5 + a_J6;

ax = a(1);
ay = a(2);
az = a(3);

dstate = [dx; dy; dz; ax; ay; az];
end