function [v1, v2] = lambert_Gauss(r1, r2, dt)
mu = 398600.44189;
TOL = 1e-8;

tm = 1;

cosdelta_theta = (dot(r1, r2)/(norm(r1)*norm(r2)));
delta_theta = asin(tm * sqrt(1 - (cosdelta_theta^2)));

l = ((norm(r1)+norm(r2)) / (4*sqrt(norm(r1)*norm(r2))*cos(delta_theta/2))) - 0.5;
m = (mu*dt^2) / ((2*sqrt(norm(r1)*norm(r2))*cos(delta_theta/2))^3);

y = 1; % initial guess
yprev = y;
ydif = 1; % random
max_its = 10000;
count = 0;

while ydif > TOL
    x = (m/(y^2)) - l;
    X = (4/3) * (1 + (6*x/5) + (48*x^2 / 35) + (480*x^3 / 315) + (5760*x^4 / 3465) + ((80640)*x^5 / (45045)) + ((1290240)*x^6 / (675675))); %+ ((23224320)*x^7 / (11486475)) + ((464486400)*x^8 / (218243025)) + ((10218700800)*x^9 / (4583103525)));
    y = 1 + X*(l + x);
    ydif = abs(y - yprev);
    yprev = y;
    count = count + 1;

    if count == max_its
        warning('Gauss method did not converge. Returned values are invalid.')
        break
    end
        
end

delta_E = asin(sqrt(x))*4;


a = ((sqrt(mu)*dt) / (2*y*sqrt(norm(r1)*norm(r2))*sin(delta_E/2)*cos(delta_theta/2)))^2;
f = 1 - ((a/norm(r1))*(1 - cos(delta_E)));
g = dt - ((a^1.5 / sqrt(mu))*(delta_E - sin(delta_E)));

p = (norm(r1)*norm(r2)*(1 - cos(delta_theta))) / (norm(r1) + norm(r2) - (2*sqrt(norm(r1)*norm(r2))*cos(delta_theta/2)*cos(delta_E/2)));
% f2 = 1 - ((norm(r2)/p)*(1-cos(delta_theta)));
% g2 = (norm(r1)*norm(r2)*sin(delta_theta)) / sqrt(mu*p);
% vtest = (r2 - (f2*r1))/g2;

v1 = (r2 - (f*r1))/g;

% fdot = sqrt(1/p)*tan(delta_theta/2)*((1-cos(delta_theta))/p - (1/norm(r2)) - (1/norm(r1)));
gdot = 1 - (norm(r1)/p)*(1 - cos(delta_theta));
v2 = ((gdot*r2)-r1)/g;
end