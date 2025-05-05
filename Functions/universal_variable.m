function [r,v] = universal_variable(r0, v0, dt, mu, TOL)

mag_r0 = sqrt(dot(r0, r0));
mag_v0 = sqrt(dot(v0, v0));
v_r0 = (dot(v0, r0))/mag_r0;
alpha = (2/mag_r0) - (mag_v0^2/mu);
univ_0 = sqrt(mu)*abs(alpha)*(dt);
z = alpha*univ_0^2;

univ = univ_0;
ratio = 1; % arbitrary
its = 0;
max_its = 10000;

while abs(ratio) > TOL && its <= max_its
its = its + 1;
C = stumpffC(z);
S = stumpffS(z);
F = (((mag_r0*v_r0)/sqrt(mu))*univ^2*C) + ((1-(alpha*mag_r0))*univ^3*S) + (mag_r0*univ - sqrt(mu)*dt);
F_prime = ((((mag_r0*v_r0)/sqrt(mu))*univ) * (1-(alpha*univ^2*S))) + ((1-alpha*mag_r0)*univ^2*C) + mag_r0;
ratio = F/F_prime;
univ = univ - ratio;
end

f = 1 - (univ^2/mag_r0)*stumpffC(alpha*univ^2);
g = dt - (1/sqrt(mu))*univ^3*stumpffS(alpha*univ^2);
r_prev = r0;
v_prev = v0;
r = f*(r_prev) + g*(v_prev);
r_mag = norm(r);
f_dot = (sqrt(mu)/(mag_r0*r_mag)) * (alpha*univ^3*stumpffS(alpha*univ^2) - univ);
g_dot = 1 - (univ^2/r_mag)*stumpffC(alpha*univ^2);
v = f_dot*(r_prev) + g_dot*(v_prev);


function c = stumpffC(z) % credit Curtis

if z > 0
c = (1 - cos(sqrt(z)))/z;
elseif z < 0
c = (cosh(sqrt(-z)) - 1)/(-z);
else
c = 1/2;
end
end

function s = stumpffS(z) % credit Curtis
if z > 0
s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0
s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
s = 1/6;
end
end

end