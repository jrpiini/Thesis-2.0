% function [sat_Pos, sat_Vel, Iterations, r1, r3] = Gauss_ODTBX(meas, time, rsite)
function [sat_Pos, sat_Vel, sat_Vel_ACO, t2, r1] = GaussExt_ODTBX(Observation, ObserverOrbit)

%% Note by Joseph Piini: this should work now, I've adapted it to take in my structures

tdif = Observation.t(end) - Observation.t(1);
t_mid = tdif/2;
[~,t2_ind] = min(abs(Observation.t-t_mid));

meas = [Observation.angles(1,:); Observation.angles(t2_ind,:); Observation.angles(end,:) ];
time = [ObserverOrbit.time_vects(1,:); ObserverOrbit.time_vects(t2_ind,:); ObserverOrbit.time_vects(end,:)];
rsite = [ObserverOrbit.r(1,:)' ObserverOrbit.r(t2_ind,:)' ObserverOrbit.r(end,:)'];

% GAUSS Gauss method of angles-only initial orbit determination
%
%    [sat_Pos, sat_Vel, Iterations] = Gauss(meas, time, rsite) Uses Gauss' 
%    initial orbit determination technique in order to provide an intertial
%    position and velocity vector for an orbiting body.
%
% INPUTS
% meas - Angle observations in a nx2 matrix [RA declination] (degrees)
% time - nx6 matrix of n datevectors for n observation times (datevector)
% rsite - 3xn matrix of site position vectors in the intertial frame (km)
%
% WARNING: Input units must be degrees for meas, datevector format
% for time, and km for rsite.
%
% OUTPUTS:
% sat_Pos - Orbiting object's position vector for second observation time (km)
% sat_Vel - Orbiting object's velocity vector for the second observation
%           time (km/s)
% Iterations -  Number of iterations to reach convergence
% r1 - Orbiting object's position vector for first observation time. Needed
%      for the initial range estimates of Double-r method. (km)
% r3 - Orbiting object's position vector for third observation time. Needed
%      for the initial range estimates of Double-r method. (km)
%
% CALLS:
% 1) LOS_Vectors
% 2) Gibbs
% 3) Herrick-Gibbs
% 4) fgcalc
% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/04/2015                Original
% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.
% disp('%%%%%%%%%% GAUSS METHOD %%%%%%%%%%');
% Gravitational Parameter
% global muglobal;
mu = 398600.44189;   % (km^3/s^2)
Ltemp = [];
for i = 1:3
    Ltemp = [Ltemp; [-cosd(meas(i,1))*cosd(meas(i,2)); sind(meas(i,1))*cosd(meas(i,2)); sind(meas(i,2))]];
    
end
LOS = LOS_Vectors(meas);
% Change in observation times
% t1 = etime(time(1,:), time(2,:));
% t3 = etime(time(3,:), time(2,:));
t1 = Observation.t(1);
t2 = Observation.t(t2_ind);
t3 = Observation.t(end);

t1 = t1-t2;
t3 = t3-t2;
% Rewriting coefficients
a1 = t3 / (t3 -t1);
a1u = t3*((t3-t1)^2 - t3^2) / (6*(t3 - t1));
a3 = -t1/ (t3-t1);
a3u = -t1*((t3-t1)^2 - t1^2) / (6*(t3 - t1));
M = LOS\ rsite;
d1 = M(2,1)*a1 - M(2,2) + M(2,3)*a3;
d2 = M(2,1)*a1u + M(2,3)*a3u;
C = dot(LOS(:,2), rsite(:,2));
% The coefficients of the eighth order polynomial that must be solved in
% order to find the middle radius magnitude can be condensed for 
% simplicity. Given the following form: r^8 + q1*r^6 + q2*r^3 + q3 = 0;
q1 = -(d1^2 + 2*C*d1 + dot(rsite(:,2),rsite(:,2)));
q2 = -2*mu*(C*d2 + d1*d2);
q3 = - mu^2*d2^2;
% The roots of the eighth order polynomial are solved using MATLAB's 
% roots() command 
% thing = [1 0 q1 0 0 q2 0 0 q3];
% eqn = @(r2) thing(1)*r2^8 + thing(2)*r2^7 + thing(3)*r2^6 + thing(4)*r2^5 + thing(5)*r2^4 + thing(6)*r2^3 + thing(7)*r2^2 + thing(8)*r2^1 + thing(9);
% eqn = @(r2) r2^8 - (r2^6*((d1^2) + (2*C*d1) + (norm(rsite(:,2))^2))) - (2*mu*r2^3*(((C*d2) + (d1*d2)))) - (mu^2 * d2^2);
% options = odeset('RelTol', 1e-8,'AbsTol',1e-8);
% root = fzero(eqn, 7300, options);
% r2 = root;
r2 = roots([1 0 q1 0 0 q2 0 0 q3]);
% Extract only the real roots
r2 = r2(r2 == real(r2));
% The positive of the two real roots represents the first iteration
% solution. Now the guess for the f and g series coefficient u can be
% updated. 
if r2(1) > 0
   r2 = r2(1);
elseif r2 (2) > 0
    r2 = r2(2);
else
    disp(mu)
    error('No positive, real roots found. No Solution.\n');
    
end
% Update guess for the f and g series ceofficient
u = mu / r2^3;
% Calculation of the original coefficients which defined the linear
% independence of each position vector
c1 = a1 + a1u*u;
c2 = -1;
c3 = a3 + a3u*u;
c = [c1; c2; c3];
% Calculate the initial guess of the slant-ranges, rho
cp = -1*(M*c);
rhonot = [cp(1)/c1; cp(2)/c2; cp(3)/c3];    % (ER)
rho = rhonot;
rho_next = [0; 0; 0];
% Initial position vectors estimate
r = [rhonot(1)*LOS(:,1) + rsite(:,1) rhonot(2)*LOS(:,2) + rsite(:,2) rhonot(3)*LOS(:,3) + rsite(:,3)];    % (ER)
% Iteration count to keep track of how many iterations until convergence
i = 0;
% Set the error tolerance for the iterative calculations of rho values
Error = rhonot(1)*1e-8;
% Counter for number of times the error tolerance is increased
error_count = 0;
%% Iterative Convergance of Solution
% Iterates until slant-range estimate stops changing significantly
while (abs(rho(1)-rho_next(1)) > Error && abs(rho(2)-rho_next(2)) > Error && abs(rho(3)-rho_next(3)) > Error)
%     disp(abs(rho(1)-rho_next(1)))
    % Check the angular separation to determine whether to use Gibbs or
    % Herrick-Gibbs method. Uses Gibbs for angles greater than 3 degrees
    % and Herrick-Gibbs for angles less than 3 degrees
    alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
    alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));
%     Use Gibbs method if separation angles are large enough
    if (abs(alpha_12) > 5 && abs(alpha_23) > 5)
        % Gibbs returns the velocity, eccentricity, and semiparameter.
        [v2, e, p] = Gibbs(r);
        
    else
        % Use Herrick-Gibbs if separation angles are small
        v2 = Herrick_Gibbs(r,time);
        % KOE = kepel(r(:,2),v2,mu);
        KOE = curtis_coe(r(:,2),v2,mu);
        p = KOE(end)*(1-KOE(2)^2);
    end
%     % Using Lambert's Problem with Minimum Energy Solution
%     v2 = Lambert_Min_Energy(r(:,2),r(:,3));
%     disp('Lambert v2')
%     disp(v2);
%     KOE = kepel(r(:,2),v2,mu);
%     p = KOE.sma*(1-KOE.ecc^2);
    
    %% Calculation of eccentricity and semiparameter using ELORB techniques
    % Angular momentum vector
    h = cross(r(:,2),v2);
    % Node vector
    n = cross([0;0;1],h);
    % Eccentricity vector
    e = ((dot(v2,v2)-mu/norm(r(:,2)))*r(:,2) - dot(r(:,2),v2)*v2)/mu;
    % Specific mechanical energy (sme)
    sme = dot(v2,v2)/2 - mu/norm(r(:,2));
    if e == 1.0
        % Semiparameter
        p = dot(h,h)/mu;
    else
        % Semimajor Axis
        a = -mu/(2*sme);
        p = a*(1-dot(e,e));
    end
   
    % Calculation of f and g coefficients
    f1 = 1 - (norm(r(:,1)) / p)*(1-cosd(-alpha_12));
    f3 = 1 - (norm(r(:,3)) / p)*(1-cosd(alpha_23));
    g1 = (norm(r(:,1))*norm(r(:,2))*sind(-alpha_12))/sqrt(mu*p);
    g3 = (norm(r(:,3))*norm(r(:,2))*sind(alpha_23))/sqrt(mu*p);
    c1 = g3 / (f1*g3 - f3*g1);
    c3 = -g1 / (f1*g3 - f3*g1);
    
    % Recalculation of slant range distances
    c = [c1; c2; c3];
    cp = M*-1*c;
    rho_temp = [cp(1)/c1; cp(2)/c2; cp(3)/c3];
   
    if i == 1
        rho = rho;
        rho_next = rho_temp;
    else
        rho = rho_next;
        rho_next = rho_temp;
    end
    i = i + 1;
%     fprintf(1,'Iteration: %f\n',i);
    r = [rho_next(1)*LOS(:,1) + rsite(:,1) rho_next(2)*LOS(:,2) + rsite(:,2) rho_next(3)*LOS(:,3) + rsite(:,3)];
    if i >= 500
       Error = Error * 10;
       error_count = error_count + 1;
       i = 0;
       r = [rhonot(1)*LOS(:,1) + rsite(:,1) rhonot(2)*LOS(:,2) + rsite(:,2) rhonot(3)*LOS(:,3) + rsite(:,3)]; 
       % fprintf(1,'Increasing error tolerance by factor of 10.\n');
    end
    if error_count == 6
       beep;
       % errordlg(['Solution diverged. Stopped after 3000 iterations'],'Computation Error','modal');
       % fprintf(1,'WARNING: Solution diverged. Stopped at 3000 iterations.\n');
       sat_Pos = [NaN NaN NaN];
       sat_Vel = [NaN NaN NaN];
       sat_Vel_ACO = [NaN NaN NaN];
       return;
    end
end
% Updated position vectors estimate
sat_Pos = r(:,2)';    % (km)
sat_Vel = v2';    % (km/s)

h_hat = cross(sat_Pos, sat_Vel) / norm(cross(sat_Pos, sat_Vel));
sat_Pos_hat = sat_Pos / norm(sat_Pos);

v2mag_cir = sqrt(mu/norm(sat_Pos));
% sat_Vel_hat = sat_Vel / norm(sat_Vel);
sat_Vel_hat = cross(h_hat, sat_Pos_hat);
sat_Vel_ACO = v2mag_cir * sat_Vel_hat;


r1 = r(:,1);
r3 = r(:,3);
Iterations = i;
end




function [LOS] = LOS_Vectors(meas)
% LOS_VECTORS Calculation of Line-Of-Site Vectors
% 
%    [LOS] = LOS_Vectors(meas) Returns the line-of-site unit vectors 
%    given topocentric right ascension and declination observations.
%
% INPUTS:
% meas - nx2 matrix containing n observations with the first column holding
%        right ascention values and the second column holding declination
%        values [RA Dec] (degrees)
%
% OUTPUTS:
% LOS - 3xn matrix of LOS vectors with each column, n, corresponding to
%       each observation given in the input
%
% CALLS:
% None
% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/03/2015                Original
% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.
LOS = [];
for i = 1:length(meas)
    RA = meas(i, 1);
    dec = meas(i, 2);
    LOS = [LOS; [cosd(dec)*cosd(RA) cosd(dec)*sind(RA) sind(dec)]];
end
LOS = LOS.';
end


function [v2, e, P] = Gibbs(r)
% GIBBS Gibbs method of velocity determination given three position vectors
%
%    [v2, e, P] = Gibbs(r) Uses the Gibbs method to return an estimation 
%    of the velocity vector for the middle observation time out of three
%    observations. The function also estimates the eccentricity, e, and
%    semilatus rectum, P.
% 
% INPUTS:
% r - 3x3 matrix containing orbiting object's position vectors at three
%     sequential observation times. (km)
%
% OUTPUTS:
% v2 - Estimation of object's velocity vector at second time (km/TU)
% e - Eccentricity estimate of the orbit
% P - Semiparameter of the orbit (km)
%
% CALLS:
% None
%
% WARNING: To use the Gibbs method, check that the position vectors are
% coplanar to a certain tolerance (<3deg) and that the separation between
% the vectors is large (>1deg). Failure to do this will cause erroneous
% estimations. This function assumes that these checks are done before
% passing arguments through, but will check again and provide error
% messages to indicate potential errors.
% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/05/2015                Original
% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.
% Gravitational Parameter
% global muglobal;
mu = 398600.44189;   % (km^3/s^2)
r1 = r(:,1);
r2 = r(:,2);
r3 = r(:,3);
Z12 = cross(r1, r2);
Z23 = cross(r2, r3);
Z31 = cross(r3, r1);
% Check that position vectors are within coplanar tolerance (~3deg)
a_COP = 90 - acosd(dot(Z23,r(:,1))/(norm(Z23)*norm(r(:,1))));
if abs(a_COP) > 3
    fprintf(1,['\nWARNING: Position Vectors are not coplanar to 3 degree', ... 
    ' tolerance!\nEstimations  using Gibbs or Herrick-Gibbs', ...
    ' may be off by a large amount.\n']); 
    fprintf(1,'Coplanar angle is: %f deg\n',a_COP); 
end 
alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));
if (abs(alpha_12) > 1 && abs(alpha_23) > 1)
else
    fprintf(1, ['\nWARNING: Angle between position vectors is smaller than',...
        ' 1 degree. This will result in erroneous results when using Gibbs',...
        ' Method!\n']);
    fprintf(1,'Angle between first and second vector is: %f deg\n',alpha_12);
    fprintf(1,'Angle between second and third vector is: %f deg\n',alpha_23);
end
N = norm(r1)*Z23 + norm(r2)*Z31 + norm(r3)*Z12;
D = Z12 + Z23 + Z31;
S = (norm(r2) - norm(r3))*r1 + (norm(r3) - norm(r1))*r2 + (norm(r1) - norm(r2))*r3;
B = cross(D, r2);
Lg = sqrt(mu/(norm(N)*norm(D)));
% Direction of eccentricity (e_hat) used to find the true anomaly
W_hat = N / norm(N);
Q_hat = S / norm(S);
e_hat = cross(Q_hat, W_hat);
v2 = (Lg/norm(r2))*B + Lg*S;    % (km/s)
e = norm(S)/norm(D);
P = norm(N)*e / norm(S);    % (km)
end


function [v2] = Herrick_Gibbs(r, time)
% HERRICK_GIBBS Herrick-Gibbs method of velocity determination given three postion vectors
% 
%    [v2] = Herrick_Gibbs(r, time)Uses the Herrick-Gibbs method to return 
%    an estimation of the velocity vector for the middle observation time 
%    out of three observations. 
%
% This method is best used with position vectors from a single pass of the
% object over a particular ground station.
%
% INPUTS:
% r - 3x3 matrix containing orbiting object's position vectors at three
%     sequential observation times. (km)
% time - 3x6 matrix of datevectors for the three observation times
%
% OUTPUTS:
% v2 - Estimation of object's velocity vector at second time (km/s)
%
% CALLS:
% None
%
% WARNING: To use the Herrick-Gibbs method, check that the position vectors are
% coplanar to a certain tolerance (<3deg) and that the separation between
% the vectors is small (< 5deg). Failure to do this will cause erroneous
% estimations. This function assumes that these checks are done before
% passing arguments through, but will check again and provide error
% messages to indicate potential errors.
% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/08/2015                Original
% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.
% Gravitational Parameter
% global muglobal;
mu = 398600.44189;   % (km^3/s^2)
r1 = r(:,1);
r2 = r(:,2);
r3 = r(:,3);
Z12 = cross(r1, r2);
Z23 = cross(r2, r3);
Z31 = cross(r3, r1);
% Check that position vectors are within coplanar tolerance (~3deg)
a_COP = 90 - acosd(dot(Z23,r(:,1))/(norm(Z23)*norm(r(:,1))));
if abs(a_COP) > 3
    fprintf(1,['\nWARNING: Position Vectors are not coplanar to 3 degree', ... 
    ' tolerance!\nEstimations  using Gibbs or Herrick-Gibbs', ...
    ' may be off by a large amount.\n']); 
    fprintf(1,'Coplanar angle is: %f deg\n',a_COP); 
end
% Check that angle between position vectors is small enough
alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));
if (abs(alpha_12) > 5 && abs(alpha_23) > 5)
    fprintf(1, ['\nWARNING: Angle between position vectors is greater than',...
        ' 3 degrees. This will result in erroneous results when using',...
        ' Herrick-Gibbs Method!\n']);
    fprintf(1,'Angle between first and second vector is: %f deg\n',alpha_12);
    fprintf(1,'Angle between second and third vector is: %f deg\n',alpha_23);
else
end
t21 = etime(time(2,:), time(1,:));
t32 = etime(time(3,:), time(2,:));
t31 = etime(time(3,:), time(1,:));
% Velocity Vector estimation using Taylor Series expansion
v2 = -t32*(1/(t21*t31)+mu/(12*norm(r(:,1))^3))*r(:,1) + ...
    (t32 - t21)*(1/(t21*t32)+mu/(12*norm(r(:,2))^3))*r(:,2)+...
    t21*(1/(t32*t31)+mu/(12*norm(r(:,3))^3))*r(:,3);    % (km/s)
end