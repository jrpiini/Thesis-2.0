%% 2LE to state vector
% Inputs two line elements, outputs state vector (similar to COEs to state)
function [r, v, h, a, theta, period] = twoline_to_state(twoLE, mu)

% Inputs
ecc = twoLE(1); % eccentricity
RAAN = twoLE(2); % right ascension of ascending node, degrees
inc = twoLE(3); % inclination, degrees
omega = twoLE(4); % argument of periapse, degrees
Me = deg2rad(twoLE(5)); % mean anomaly, degrees
n = twoLE(6); % mean motion, revs/day

n = (n/(24*60*60)*(2*pi)); % revs/day --> rads/sec
a = (mu/(n^2))^(1/3); % semi-major axis
period = (2*pi) * sqrt((a^3)/mu);
h = sqrt(mu*a*(1-(ecc)^2)); % specific angular momentum

E_eqn = @(E) E - ecc*sin(E) - Me; % Eccentric anomaly
[bracket_matrix] = bracketing_algorithm(E_eqn, -50, .1, 1);
x_0 = bracket_matrix(1,1);
x_1 = bracket_matrix(1,2);
[z, ~, ~] = bisection_algorithm(x_0, x_1, E_eqn, 1e-10);
E = z;
theta = 2*atand((tan(E/2) * sqrt((1+ecc)/(1-ecc)))); % true anomaly

if theta > 360 % adjust if theta is not between 0 and 360
    multiples = floor(theta/360);
    theta = theta - (360*multiples);

    elseif theta < 0
    multiples = ceil(abs(theta)/360);
    theta = theta + (360*multiples);
end

% perifocal
rp = (h^2/mu) * (1/(1 + ecc*cosd(theta))) * (cosd(theta)*[1;0;0] + sind(theta)*[0;1;0]);
vp = (mu/h) * (-sind(theta)*[1;0;0] + (ecc + cosd(theta))*[0;1;0]);

R3_W = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
R1_i = [1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
R3_w = [cosd(omega) sind(omega) 0; -sind(omega) cosd(omega) 0; 0 0 1];

Q = (R3_w*R1_i*R3_W)';

% Final vectors
r = (Q*rp)';
v = (Q*vp)';





%% Bisection and bracketing root-finding functions

% Bisection method; root-finder used to find eccentric anomaly
function [r, count, output] = bisection_algorithm(a, b, p, TOL)
output = [a, b, p((a+b)/2) (a+b)/2 (b-a)/2]; 
count = 1;
max_its = 10000; 
while ((b - a)/2 > TOL)
    count = count + 1; 
    if count == max_its
        break
    end 
    c = (a + b)/2; 
    if p(c) == 0   
        break;
    end    
    if (p(a)*p(c) < 0) 
        b = c; 
    else
        a = c; 
    end
    output = [output; [a b p((a+b)/2) (a+b)/2 (b-a)/2]]; 
end
r = (a+b)/2;
end

end

% Bracketing algorithm
function [matrix] = bracketing_algorithm(p, a, h, r_expected)
initial_guess = a;
matrix = zeros(r_expected, 2);
k = 0;
n = 0;
max_its = 10000;
b = a + h;

while n < r_expected 
    f = p(a);
    while f*p(b) > 0 
    b = b + h;
    k = k + 1;
        if k == max_its
        break
        end        
    end   
    if p(a)*p(b) < 0   
        n = n+1; 
        a = b - h;
        matrix(n,:) = [a, b];
        a = b; 
        b = a+h; 
    end
     
end

if n < r_expected
    a = initial_guess;
    b = a - h;
    f = p(a);
    while f*p(b) > 0
    b = b - h;   
    k = k + 1;
        if k == max_its
        break
        end       
    end   
    if p(a)*p(b) < 0     
        n = n+1; 
        a = b + h;
        matrix(n,:) = [a, b]; % return a and b values in matrix
        a = b;
        b = a-h;
    end
     
end


end



