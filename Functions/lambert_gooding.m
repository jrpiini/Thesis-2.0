function [vi, vf] = lambert_gooding(r1, r2, dt)
mu = 398600.44189;
ot = 1;
% this section contains the functions for the gooding lambert solution
% inputs
%  mu = gravitational constant (kilometers^3/seconds^2)
%  r1 = initial state vector (kilometers)
%  r2 = final state vector (kilometers)
%  dt = transfer time (seconds)
%  ot = orbit type (1 = prograde, 2 = retrograde)
% outputs
%  v1 = initial velocity vector of transfer orbit (kilometers/second)
%  v2 = final velocity vector of transfer orbit (kilometers/second)
%  R. H. Gooding, Technical Report 88027
%  On the Solution of Lambert's Orbital Boundary-Value Problem,
%  Royal Aerospace Establishment, April 1988
%  R. H. Gooding, Technical Memo SPACE 378
%  A Procedure for the Solution of Lambert's Orbital Boundary-Value Problem
%  Royal Aerospace Establishment, March 1991
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize time difference
tdiff = 1.0e8;
nrev = 0;
itermax = 20;
tol = 1.0e-10;
error1 = 0;
r1mag = norm(r1);
r2mag = norm(r2);

% if ~isreal(r1mag) || ~isreal(r2mag)
%     vi = [0 0 0]; vf = vi;
%     return
% end

r1unit = r1 / r1mag;
r2unit = r2 / r2mag;
c12 = cross(r1, r2);
c12mag = norm(c12);
c12unit = c12 / c12mag;
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
% tangential unit vectors
tan1u = cross(c12unit, r1unit);
tan2u = cross(c12unit, r2unit);
if (nu >= pi)
    tan1u = -tan1u;
    tan2u = -tan2u;
end
% constants determined from problem geometry
c = sqrt(r1mag * r1mag + r2mag * r2mag - 2.0 * r1mag * r2mag * cos(nu));
s = (r1mag + r2mag + c) / 2.0;
t = sqrt(8.0 * mu / s^3) * dt;
q = sqrt(r1mag * r2mag) / s * cos(nu * 0.5);
if (nu <= pi)
    q = abs(q);
else
    q = -abs(q);
end
% initial conditions
[t0, ~, ~] = time_derivatives(0.0, q, nrev);
if ((t0 - t) < 0.0)
    signx = -1;
else
    signx = 1;
end
% constants for initial conditions splicing
c1 = 0.5;
c2 = 0.03;
% initial values for single revolution case
if (signx == 1)
    x0 = t0 * (t0 - t) / (4.0 * t);
else
    x01 = -(t - t0) / (t - t0 + 4.0);
    x02 = -sqrt((t - t0) / (t + 0.5 * t0));
    w_big = x01 + c1 * sqrt(2.0 - nu / pi);
    if (w_big >= 0.0)
        x03 = x01;
    else
        w_little = sqrt(sqrt(sqrt(sqrt(-w_big))));
        x03 = x01 + w_little * (x02 - x01);
    end
    lam = 1.0 + c1 * x03 * (4.0 / (4.0 - t - t0)) - c2 * x03^2 ...
        * sqrt(1.0 + x01);
    x0 = lam * x03;
end
if (x0 <= -1)
    error1 = 1;
end
% solve each case given the initial conditions above for 0-revolution case
gamma = sqrt(mu * s / 2.0);
if (error1 == 1)
    % no solution found
    v1(1:3) = 1.0e7;
    v2(1:3) = 1.0e7;
else
    % initialize iteration
    iter = 0;
    x = x0;
    while (abs(tdiff) >= tol && iter <= itermax)
        iter = iter + 1;
        if ~isreal(x)
            vi = [0 0 0]; vf = vi;
            return
        end
        [t2, td1, td2] = time_derivatives(x, q, nrev);
        tdiff = t2 - t;
        x_new = x - 2.0 * tdiff * td1 / (2.0 * td1 * td1 - tdiff * td2);
        x = x_new;
    end
    if (iter >= itermax)
        % no solution found
        v1(1:3) = 1.0e12;
        v2(1:3) = 1.0e12;
    else
        if (c == 0)
            sig = 1.0;
            rho = 0.0;
            z = abs(x);
        else
            sig = 2.0 * sqrt(r1mag * r2mag / c^2) * sin(0.5 * nu);
            rho = (r1mag - r2mag) / c;
            z = sqrt(1.0 + q^2 * (x^2 - 1.0));
        end
        % radial components of velocity (kilometers/second)
        dum = gamma * ((q * z - x) - rho * (q * z + x)) / r1mag;
        vr1 = dum * r1unit;
        dum = -gamma * ((q * z - x) + rho * (q * z + x)) / r2mag;
        vr2 = dum * r2unit;
        % tangential components of velocity (kilometers/second)
        dum = gamma * sig * (z + q * x);
        vt1 = dum / r1mag * tan1u;
        vt2 = dum / r2mag * tan2u;
        v1 = vr1 + vt1;
        v2 = vr2 + vt2;
    end
end
% load solution velocity vectors
vi = v1;
vf = v2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, td1, td2] = time_derivatives(x, q, nrev)
% this function calculates the time-of-flight equations and the derivatives
% inputs
%  x    = current x estimate
%  q    = lambert parameter
%  nrev = number of complete revolutions
% outputs
%  t   = calculated time-of-flight
%  td1 = time-of-flight first derivative
%  td2 = time-of-flight second derivative
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1.0e-8;
delt = 1.0e-1;
xx = x;
if (xx < -1.0)
    xx = abs(x) - 2.0;
elseif (x == -1)
    xx = x + delt;
end
e = xx * xx - 1.0;
if (xx == 1)
    % exact parabolic solution, probably never actually used
    t = 4.0 / 3.0 * (1.0 - q^3);
    td1 = 4.0 / 5.0 * (q^5 - 1.0);
    td2 = td1 + 120.0 / 70.0 * (1.0 - q^7);
elseif (abs(xx - 1.0) < tol)
    % near parabolic solution, use the trans/series representation
    [sig1, dsig1, d2sig1] = sigma(-e);
    [sig2, dsig2, d2sig2] = sigma(-e * q * q );
    t = sig1 - q^3 * sig2;
    td1 = 2.0 * x * (q^5.0 * dsig2 - dsig1);
    td2 = td1 / x + 4.0 * x * x * (d2sig1 - q^7.0 * d2sig2);
else
    % all cases not exactly parablolic or close to parabolic
    y = sqrt(abs(e));
    z = sqrt(1.0 - q^2 + q^2 * x^2);
    f = y * (z - q * x);
    u = -e;
    % beta = q * z - x;
    g = (x^2 - q^2 * (-e)) / (x * z - q * (-e));

    % if ~isreal(f) || ~isreal(g)
    %     asdftg = 5;
    % end

    if (e < 0)
        d = atan2(f, g) + pi * nrev;
    else
        d = atanh(f / g);
        % d = log(f+g)
    end
    t = 2.0 * (x - q * z - d / y) / e;
    td1 = (3.0 * x * t + 4.0 * q^3 * x / z - 4.0) / u;
    td2 = (3.0 * t + 5.0 * x * td1 + 4.0 * (q / z)^3 * (1.0 - q^2)) / u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigm, dsigma, d2sigma] = sigma(u)
% this function computes the sigma and sigma derivatives for the
% gooding lambert solution
% input
%  u = initial u value for the sigma function
% output
%  sigm    = sigma(u)
%  dsigma  = d(sigma)/du
%  d2sigma = d2(sigma)/du2
% Orbital Mechanics with MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigm = -2.0 * (sqrt(u) * sqrt(1.0 - u) - asin(sqrt(u))) / sqrt(u)^3;
dsigma = -(u / sqrt(1.0 - u) + (3.0 * asin(sqrt(u))) / sqrt(u) ...
    - 3.0 / sqrt((1.0 - u))) / u^2;
d2sigma = (15.0 * asin(sqrt(u)) * sqrt(1.0 - u)^3 - 15.0 * sqrt(u) ...
    + 20.0 * sqrt(u)^3 - 3.0 * sqrt(u)^5) / (2.0 * sqrt(u)^7 ...
    * sqrt((1.0 - u))^3);
