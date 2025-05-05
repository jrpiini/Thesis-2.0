function [vi, vf] = lambert_battin(r1, r2, dt)
mu = 398600.44189;
ot = 1;
% this function contains the battin lambert solution method
% inputs
%  mu = gravitational constant (kilometers^3/seconds^2)
%  r1 = initial position vector (kilometers)
%  r2 = final position vector (kilometers)
%  dt = transfer time (seconds)
%  ot = orbit type (1 = prograde, 2 = retrograde)
% outputs
% vi = initial velocity vector of transfer orbit (kilometers/second)
% vf = final velocity vector of transfer orbit (kilometers/second)
% Battin, R. An Introduction to the Mathematics and Methods of Astrodynamics,
% Chapter 7: Solving Lambert's Problem, AIAA Education Series,
% 1801 Alexander Bell Drive, Reston, VA, Revised Edition edition, 1999
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convergence tolerance
tol = 1.0e-8;
r1mag = norm(r1);
r2mag = norm(r2);
% determine true anomaly angle here (radians)
c12 = cross(r1, r2);
nu = acos(dot(r1, r2) / (r1mag * r2mag));
% determine the true anomaly angle using the orbit type
% 1 is prograde, 2 is retrograde
if (ot == 1)
    if (c12(3) <= 0.0)
        nu = 2.0 * pi - nu;
    end
end
if (ot == 2)
    if (c12(3) >= 0.0)
        nu = 2.0 * pi - nu;
    end
end
c = sqrt(r1mag * r1mag + r2mag * r2mag - 2.0 * r1mag * r2mag * cos(nu));
s = (r1mag + r2mag + c) / 2.0;
eps = (r2mag - r1mag) / r1mag;
lam = sqrt(r1mag * r2mag) * cos(nu * 0.5) / s;
t = sqrt(8.0 * mu / s^3) * dt;
t_p = 4.0 / 3.0 * (1.0 - lam^3);
m = t^2 / (1.0 + lam)^6;
tansq2w = (eps * eps * 0.25) / (sqrt(r2mag / r1mag) + r2mag / r1mag ...
    * (2.0 + sqrt(r2mag / r1mag)));
rop = sqrt(r2mag * r1mag) * (cos(nu * 0.25) * cos(nu * 0.25) + tansq2w);
if (nu < pi)
    ltop = (sin(nu * 25.0e-2) * sin(nu * 25.0e-2) + tansq2w);
    l = ltop / (ltop + cos(nu * 5.0e-1));
else
    ltop = cos(nu * 25.0e-2) * cos(nu * 25.0e-2) + tansq2w;
    l = (ltop - cos(nu * 5.0e-1)) / ltop;
end
% initial guess is set here
if (t <= t_p)
    x = 0.0;
else
    x = l;
end
dx = 1.0;
iter = 0;
itermax = 20;
% this loop does the successive substitution
while (dx >= tol && iter <= itermax)
    xi = xi_battin(x);
    denom = (1.0 + 2.0 * x + l) * (4.0 * x + xi * (3.0 + x));
    h1 = (l + x)^2 * (1.0 + 3.0 * x + xi) / denom;
    h2 = (m * (x - l + xi)) / denom;
    b = 27.0 * h2 * 25.0e-2 / (1.0 + h1)^3;
    u = b / (2.0 * (sqrt(1.0 + b) + 1.0));
    k = k_battin(u);
    y = (1.0 + h1) / 3.0 * (2.0 + sqrt(1.0 + b) / (1.0 + 2.0 * u * k * k));
    xnew = sqrt(((1.0 - l) / 2.0)^2 + m / (y * y)) - (1.0 + l) / 2.0;
    dx = abs(x - xnew);
    x = xnew;
    iter = iter + 1;
end
if (iter >= itermax)
    % if a solution isn't found the final velocities are
    % output as a large number
    v1(1:3) = 1.0e12;
    v2(1:3) = 1.0e12;
else
    a = mu * dt * dt / (16.0 * rop * rop * x * y * y);
    fg = fg_battin(mu, a, s, c, nu, dt, r1mag, r2mag);
    v1 = (r2 - fg(1) * r1) / fg(2);
    v2 = (fg(3) * r2 - r1) / fg(2);
end
% load velocity vectors solution
vi = v1;
vf = v2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function xi = xi_battin(x)
% this function computes the first continued fraction for the battin
% algorithm
% input
%  x = current x estimate
% output
%  xi = value for the continued fraction
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = [0.25396825396825395, 0.25252525252525254, 0.25174825174825177, ...
    0.25128205128205128, 0.25098039215686274, 0.25077399380804954, ...
    0.25062656641604009, 0.25051759834368531, 0.25043478260869567, ...
    0.25037037037037035, 0.25031928480204341, 0.25027808676307006, ...
    0.25024437927663734, 0.25021645021645023, 0.25019305019305021, ...
    0.25017325017325015, 0.25015634771732331, 0.25014180374361883, ...
    0.25012919896640828, 0.25011820330969264];
eta = x / (sqrt(1.0 + x) + 1.0)^2;
xi = 8.0 * (sqrt(1.0 + x) + 1.0) / (3.0 + 1.0 / ...
    (5.0 + eta + 9.0 / 7.0 * eta / (1.0 + c(1) * eta / ...
    (1.0 + c(2) * eta / (1.0 + c(3) * eta / (1.0 + c(4) * eta / ...
    (1.0 + c(5) * eta / (1.0 + c(6) * eta / (1.0 + c(7) * eta / ...
    (1.0 + c(8) * eta /(1.0 + c(9) * eta / (1.0 + c(10) * eta / ...
    (1.0 + c(11) * eta / (1.0 + c(12) * eta / (1.0 + c(13) * eta /...
    (1.0 + c(14) * eta / (1.0 + c(15) * eta / (1.0 + c(16) * eta / ...
    (1.0 + c(17) * eta / (1.0 + c(18) * eta / (1.0 + c(19) * eta / ...
    (1.0 + c(20) * eta))))))))))))))))))))));
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
function k = k_battin(u)
% this function computed the k continued fraction for the battin
% lambert solution algorithm
% input
%  u = same as the u variable in the battin description
% output
%  k = continued fraction value
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = [0.33333333333333331, 0.14814814814814814, 0.29629629629629628, ...
    0.22222222222222221, 0.27160493827160492, 0.23344556677890010, ...
    0.26418026418026419, 0.23817663817663817, 0.26056644880174290, ...
    0.24079807361541108, 0.25842383737120578, 0.24246606855302508, ...
    0.25700483091787441, 0.24362139917695474, 0.25599545906059318, ...
    0.24446916326782844, 0.25524057782122300, 0.24511784511784512, ...
    0.25465465465465464, 0.24563024563024563, 0.25418664443054689];
k = d(1) / (1.0 + d(2) * u / (1.0 + d(3) * u / (1.0 + d(4) * u /...
    (1.0 + d(5) * u / (1.0 + d(6) * u/ (1.0 + d(7) * u / ...
    (1.0 + d(8) * u / (1.0 + d(9) * u / (1.0 + d(10) * u / ...
    (1.0 + d(11) * u / (1.0 + d(12) * u / (1.0 + d(13) * u / ...
    (1.0 + d(14) * u / (1.0 + d(15) * u / (1.0 + d(16) * u / ...
    (1.0 + d(17) * u / (1.0 + d(18) * u / (1.0 + d(19) * u / ...
    (1.0 + d(20) * u / (1.0 + d(21) * u))))))))))))))))))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fg = fg_battin(mu, a, s, c, nu, t, r1, r2)
% this function computes the lagrange coefficient values to compute
% the final 2 velocity vectors
% input
%  mu = gravitational constant (kilometers^3/seconds^2)
%  a  =
%  s  =  semiparameter
%  c  =  chord
%  nu =  true anomaly angle (radians),
%  t  =  scaled time-of-flight parameter
%  r1 =  initial radius magnitude (kilometers)
%  r2 =  final radius magnitude (kilometers)
% output
%  fg(1) = f lagrange coefficient
%  fg(2) = g lagrange coefficient
%  fg(3) = gdot lagrange coefficient
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg = zeros(3, 1);
small_number = 1.0e-3;
if (a > small_number)
    be = 2.0 * asin((sqrt((s - c) / (2.0 * a))));
    if (nu > pi)
        be = -be;
    end
    a_min = s * 5.0e-1;
    t_min = sqrt(a_min^3 / mu) * (pi - be + sin(be));
    dum = (sqrt(s / (2.0 * a)));
    ae = 2.0 * asin(dum);
    if (t > t_min)
        ae = 2.0 * pi - ae;
    end
    de = ae - be;
    f = 1.0 - a / r1 * (1.0 - cos(de));
    g = t - sqrt(a * a * a / mu) * (de - sin(de));
    gdot = 1.0 - a / r2 * (1.0 - cos(de));
elseif (a < -small_number)
    ah = 2.0 * asinh(sqrt(s / (-2.0 * a)));
    bh = 2.0 * asinh(sqrt(( s - c) / (-2.0 * a)));
    dh = ah - bh;
    f = 1.0 - a / r1 * (1.0 - cosh(dh));
    g = t - sqrt(-a^3 / mu) * (sinh(dh) - dh);
    gdot = 1.0 - a / r2 * (1.0 - cosh(dh));
else
    f = 0.0;
    g = 0.0;
    gdot = 0.0;
end
% load lagrange coefficients
fg(1) = f;
fg(2) = g;
fg(3) = gdot;
