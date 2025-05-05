%% Testing out using a non-random initialization
% clear
% close all
% clc
delete(gcp('nocreate'))

%% Choose MC Options

% noise = 0; % arcseconds (desired noise to add to angles)
num_its = 300; % MC iterations
spacing = NaN; % ~3 Hz observation frequency
n_obs = 15; % fixing at 15 for varying arc lengths
% arc_length = 120; % desired arc length, seconds ! must be less than end time

Lambert_Type = 1; % 1 = Gauss, 2 = UV, 3 = Gooding, 4 = IzzoGooding, 5 = Battin


%%
% virtual = input('Using virtual computer labs? (1 for yes, 0 for no): ');
% ifVirtual = logical(virtual);
ifVirtual = boolean(0);


mu = 398600.44189;
TOL = 1e-8;

options = odeset('AbsTol',1e-8, 'RelTol',1e-8);

%% Making fake data

load time_span_6_58.mat % actual JD of observations
right_ascension = 0.*time_span; declination = 0.*time_span; % for sake of calling function
time_span = time_span(1:end)';

if ~isnan(spacing)
% if fixing frequency
time_span = (86400*time_span(1)):spacing:(86400*time_span(1) + arc_length);
time_span = (time_span ./ 86400)';
else
% if fixing # of obs
time_span = (linspace(time_span(1), time_span(1)+(arc_length/86400), n_obs))';
end

% Object 5730
% 1  5730U 71119B   18320.42753668 +.00001327 +00000-0 +88563-4 0  9996
% 2  5730 073.8926 302.5269 0692673 081.9139 285.9934 14.03000486278219
ds = '18320.42753668';
yeardayhour = str2double(regexp(ds, '(\d{2})(\d{3})(\.\d+)', 'tokens', 'once'));
dn = datenum(yeardayhour(1) + 2000, 0, yeardayhour(2), 24 * yeardayhour(3), 0, 0);
dt = datetime(dn, 'ConvertFrom', 'datenum');
TLE_JD = juliandate(dt);
[r_true, v_true] = twoline_to_state([.0692673, 302.5269, 73.8926, 81.9139, 285.9934, 14.03000486278219], mu); % at TLE time
[tout, statenew] = ode45(@cowell, [0 86400*(time_span(1)-TLE_JD)], [r_true v_true], options);
TargetOrbit.r0 = statenew(end,1:3); TargetOrbit.v0 = statenew(end,4:6);

TargetOrbit.r(1,:) = TargetOrbit.r0;
TargetOrbit.v(1,:) = TargetOrbit.v0;
TargetOrbit.t(1) = 0;
ObserverOrbit.t(1) = 0;

mu = 398600.44189;
phi_gd = 37.1384; % lat
lambda = -122.2110; % long
f = (6378-6357)/6378;
h_ellp = 684/1000; % alt, km

[~, statenew_target] = ode45(@cowell, 86400*(time_span-time_span(1)), [TargetOrbit.r0, TargetOrbit.v0], options);
TargetOrbit.r = statenew_target(:,1:3);
TargetOrbit.v = statenew_target(:,4:6);
for i = 1:height(time_span)
    % [TargetOrbit.r(i,:), TargetOrbit.v(i,:)] = universal_variable(TargetOrbit.r0, TargetOrbit.v0, (time_span(i)-time_span(1))*86400, mu, TOL);
    ObserverOrbit.t(i,:) = (time_span(i)-time_span(1))*86400;
    r_S_ECI = lla2eci([phi_gd lambda h_ellp*1000], [(year(datetime(time_span(i), 'ConvertFrom', 'juliandate'))), (month(datetime(time_span(i), 'ConvertFrom', 'juliandate'))), (day(datetime(time_span(i), 'ConvertFrom', 'juliandate'))), (hour(datetime(time_span(i), 'ConvertFrom', 'juliandate'))), (minute(datetime(time_span(i), 'ConvertFrom', 'juliandate'))), (second(datetime(time_span(i), 'ConvertFrom', 'juliandate')))]); % m
    ObserverOrbit.r(i,:) = r_S_ECI./1000; % km
    Observation.r_site(i,:) = r_S_ECI./1000; % km
    
    ObserverOrbit.time_vects(i,:) = JD_to_UTC(time_span(i)); % needed for Gauss extended
end

for i = 1:height(Observation.r_site)
    rho = TargetOrbit.r(i,:) - Observation.r_site(i,:);
    Observation.angles(i,:) = rho2RaDec_topo(rho);
    Observation.angles(i,1) = Observation.angles(i,1);% + (randn*noise/3600); % adding 1 arcsecond error noise to measurements
    Observation.angles(i,2) = Observation.angles(i,2);% + (randn*noise/3600); % adding 1 arcsecond error noise to measurements
    Observation.LOS_measurements(i,:) = LOS_from_RADec(Observation.angles(i,1), Observation.angles(i,2));
end
Observation.t = ObserverOrbit.t;

tdif = Observation.t(end) - Observation.t(1);
t_mid = tdif/2;
[~,t2_ind] = min(abs(Observation.t-t_mid));

%% Begin Monte Carlo

% RA0 = Observation.angles(1,1); Dec0 = Observation.angles(1,2);
% RAf = Observation.angles(end,1); Decf = Observation.angles(end,2);
Observation0 = Observation;

fileName = ['MonteCarlo_' + string(year(now)) + '_' + string(month(now)) + '_' + string(day(now))+ '_' + string(hour(now)) + '_' + string(minute(now)) + '.mat'];
if ~ifVirtual
fileDir = 'C:\Users\josep\OneDrive - Cal Poly\Desktop\Thesis Work\Test Cases\Made Up Data\LEO\Highish Ecc\Nonrandom Lambert Monte Carlo Results\Made Up (Updated)\New ACO';
else
% fileDir = 'C:\ProgramData\UserDataFolders\S-1-5-21-3066599673-3616497067-1000726541-1012\My Files\Temporary Files'; % If using virtual computer labs
fileDir = 'C:\ProgramData\UserDataFolders\S-1-5-21-763851881-3117066521-3573464578-1012\My Files\OneDrive\Files\Desktop\Thesis Work\Test Cases\Made Up Data\LEO\Highish Ecc\Lambert Monte Carlo Results'; % If using virtual
end
full = fullfile(fileDir, fileName);


hours = 12;
parpool('IdleTimeout', hours*60); % idle time of 12 hours

rng("shuffle")
% added_noise = randn(size(Observation0.angles));

tic
parfor(i = 1:num_its, Inf)
% i = 1;
    % Perturb all measurements
    Observation(i) = Observation0;
    Observation(i).angles = Observation(i).angles + ((noise/3600)*randn(size(Observation(i).angles)));
    % Observation(i).angles = Observation(i).angles + ((noise/3600)*added_noise);
    Observation(i).LOS_measurements = LOS_from_RADec(Observation(i).angles(:,1), Observation(i).angles(:,2));

    % Genetic Algorithm (all obs points)
    % tic
    [chromosomes(i), Elite(i)] = GA_IOD_lambert_nonrandom(Observation(i), Lambert_Type, boolean(0), boolean(0));
    MonteCarlo(i).r0 = Elite(i).r0(1,:);
    MonteCarlo(i).v0 = Elite(i).v0(1,:);
    MonteCarlo(i).ranges = Elite(i).ranges(1,:);
    MonteCarlo(i).run_time = chromosomes(i).run_time(1,:);
    MonteCarlo(i).coes_t0 = curtis_coe(Elite(i).r0(1,:), Elite(i).v0(1,:), mu);
    % toc

    % ACO Genetic Algorithm (all obs points)
    [chromosomes_ACO(i), Elite_ACO(i)] = GA_IOD_lambert_nonrandom(Observation(i), Lambert_Type, boolean(0), boolean(1));
    MonteCarlo(i).r0_ACO = Elite_ACO(i).r0(1,:);
    MonteCarlo(i).v0_ACO = Elite_ACO(i).v0(1,:);
    MonteCarlo(i).ranges_ACO = Elite_ACO(i).ranges(1,:);
    MonteCarlo(i).run_time_ACO = chromosomes_ACO(i).run_time(1,:);
    MonteCarlo(i).coes_t0_ACO = curtis_coe(Elite_ACO(i).r0(1,:), Elite_ACO(i).v0(1,:), mu);

    % Genetic Algorithm (only 3 points)
    % tic
    [chromosomes_3pts(i), Elite_3pts(i)] = GA_IOD_lambert_nonrandom(Observation(i), Lambert_Type, boolean(1), boolean(0));
    MonteCarlo(i).r0_3pts = Elite_3pts(i).r0(1,:);
    MonteCarlo(i).v0_3pts = Elite_3pts(i).v0(1,:);
    MonteCarlo(i).ranges_3pts = Elite_3pts(i).ranges(1,:);
    MonteCarlo(i).run_time_3pts = chromosomes_3pts(i).run_time(1,:);
    MonteCarlo(i).coes_t0_3pts = curtis_coe(Elite_3pts(i).r0(1,:), Elite_3pts(i).v0(1,:), mu);
    % toc

    % % ACO Genetic Algorithm (only 3 points)
    % [chromosomes_3pts_ACO(i), Elite_3pts_ACO(i)] = GA_IOD_lambert_nonrandom(Observation(i), Lambert_Type, boolean(1), boolean(1));
    % MonteCarlo(i).r0_3pts_ACO = Elite_3pts_ACO(i).r0(1,:);
    % MonteCarlo(i).v0_3pts_ACO = Elite_3pts_ACO(i).v0(1,:);
    % MonteCarlo(i).ranges_3pts_ACO = Elite_3pts_ACO(i).ranges(1,:);
    % MonteCarlo(i).run_time_3pts_ACO = chromosomes_3pts_ACO(i).run_time(1,:);
    % MonteCarlo(i).coes_t0_3pts_ACO = curtis_coe(Elite_3pts_ACO(i).r0(1,:), Elite_3pts_ACO(i).v0(1,:), mu);
   
    % Gauss extended / ACO
    % tic
    [gauss_ext(i).r2, gauss_ext(i).v2, gauss_ext_ACO(i).v2, gauss_ext(i).t2] = GaussExt_ODTBX(Observation(i), ObserverOrbit);
    % toc

    % % Gauss extended ACO
    % [gauss_ext_ACO(i).r2, gauss_ext_ACO(i).v2, gauss_ext_ACO(i).t2] = GaussExtACO_ODTBX(Observation(i), ObserverOrbit);

    % % Double-r
    % [doubleR(i).r2, doubleR(i).v2, doubleR(i).t2] = AERO557doubleR(Observation(i), ObserverOrbit);

end
t_MC = toc;

%% Reformat structure for convenience
for i = 1:length(MonteCarlo)
    % GA (all obs points)
    MonteCarlo_mod.r0_GA(i,:) = MonteCarlo(i).r0;
    MonteCarlo_mod.v0_GA(i,:) = MonteCarlo(i).v0;
    MonteCarlo_mod.ranges_GA(i,:) = MonteCarlo(i).ranges;
    MonteCarlo_mod.a_t0_GA(i,:) = MonteCarlo(i).coes_t0(end);
    MonteCarlo_mod.ecc_t0_GA(i,:) = MonteCarlo(i).coes_t0(2);
    MonteCarlo_mod.inc_t0_GA(i,:) = MonteCarlo(i).coes_t0(4);
    MonteCarlo_mod.RAAN_t0_GA(i,:) = MonteCarlo(i).coes_t0(3);
    MonteCarlo_mod.AoP_t0_GA(i,:) = MonteCarlo(i).coes_t0(5);
    MonteCarlo_mod.TA_t0_GA(i,:) = MonteCarlo(i).coes_t0(6);

    % GA (3 obs points)
    MonteCarlo_mod.r0_GA_3pts(i,:) = MonteCarlo(i).r0_3pts;
    MonteCarlo_mod.v0_GA_3pts(i,:) = MonteCarlo(i).v0_3pts;
    MonteCarlo_mod.ranges_GA_3pts(i,:) = MonteCarlo(i).ranges_3pts;
    MonteCarlo_mod.a_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(end);
    MonteCarlo_mod.ecc_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(2);
    MonteCarlo_mod.inc_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(4);
    MonteCarlo_mod.RAAN_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(3);
    MonteCarlo_mod.AoP_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(5);
    MonteCarlo_mod.TA_t0_GA_3pts(i,:) = MonteCarlo(i).coes_t0_3pts(6);

    % ACO GA (all obs points)
    % h_hat = (cross(MonteCarlo(i).r0, MonteCarlo(i).v0) / norm(cross(MonteCarlo(i).r0, MonteCarlo(i).v0))); r_hat = MonteCarlo(i).r0/norm(MonteCarlo(i).r0); v2mag_cir = sqrt(mu/norm(MonteCarlo(i).r0)); v_ACO_hat = cross(h_hat, r_hat); v_ACO = v2mag_cir * v_ACO_hat;
    % MonteCarlo_mod.r0_GA_ACO(i,:) = MonteCarlo(i).r0;
    % MonteCarlo_mod.v0_GA_ACO(i,:) = v_ACO;
    % coe = curtis_coe(MonteCarlo_mod.r0_GA_ACO(i,:), MonteCarlo_mod.v0_GA_ACO(i,:), mu);
    % % MonteCarlo_mod.ranges_GA_ACO(i,:) = MonteCarlo(i).ranges_ACO;
    % MonteCarlo_mod.a_t0_GA_ACO(i,:) = coe(end);
    % MonteCarlo_mod.ecc_t0_GA_ACO(i,:) = coe(2);
    % MonteCarlo_mod.inc_t0_GA_ACO(i,:) = coe(4);
    % MonteCarlo_mod.RAAN_t0_GA_ACO(i,:) = coe(3);
    % MonteCarlo_mod.AoP_t0_GA_ACO(i,:) = coe(5);
    % MonteCarlo_mod.TA_t0_GA_ACO(i,:) = coe(6);
    MonteCarlo_mod.r0_GA_ACO(i,:) = MonteCarlo(i).r0_ACO;
    MonteCarlo_mod.v0_GA_ACO(i,:) = MonteCarlo(i).v0_ACO;
    % MonteCarlo_mod.ranges_GA(i,:) = MonteCarlo(i).ranges;
    MonteCarlo_mod.a_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(end);
    MonteCarlo_mod.ecc_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(2);
    MonteCarlo_mod.inc_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(4);
    MonteCarlo_mod.RAAN_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(3);
    MonteCarlo_mod.AoP_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(5);
    MonteCarlo_mod.TA_t0_GA_ACO(i,:) = MonteCarlo(i).coes_t0_ACO(6);

    % ACO GA (3 obs points)
    h_hat = (cross(MonteCarlo(i).r0_3pts, MonteCarlo(i).v0_3pts) / norm(cross(MonteCarlo(i).r0_3pts, MonteCarlo(i).v0_3pts))); r_hat = MonteCarlo(i).r0_3pts/norm(MonteCarlo(i).r0_3pts); v2mag_cir = sqrt(mu/norm(MonteCarlo(i).r0_3pts)); v_ACO_hat = cross(h_hat, r_hat); v_ACO = v2mag_cir * v_ACO_hat;
    MonteCarlo_mod.r0_GA_3pts_ACO(i,:) = MonteCarlo(i).r0_3pts;
    MonteCarlo_mod.v0_GA_3pts_ACO(i,:) = v_ACO;
    coe = curtis_coe(MonteCarlo_mod.r0_GA_ACO(i,:), MonteCarlo_mod.v0_GA_ACO(i,:), mu);
    % MonteCarlo_mod.ranges_GA_3pts_ACO(i,:) = MonteCarlo(i).ranges_3pts_ACO;
    MonteCarlo_mod.a_t0_GA_3pts_ACO(i,:) = coe(end);
    MonteCarlo_mod.ecc_t0_GA_3pts_ACO(i,:) = coe(2);
    MonteCarlo_mod.inc_t0_GA_3pts_ACO(i,:) = coe(4);
    MonteCarlo_mod.RAAN_t0_GA_3pts_ACO(i,:) = coe(3);
    MonteCarlo_mod.AoP_t0_GA_3pts_ACO(i,:) = coe(5);
    MonteCarlo_mod.TA_t0_GA_3pts_ACO(i,:) = coe(6);

    % Gauss Extended
    [gauss_ext(i).r0, gauss_ext(i).v0] = universal_variable(gauss_ext(i).r2, gauss_ext(i).v2, (ObserverOrbit.t(1)-gauss_ext(i).t2), mu, TOL);
    % [~, statenew_GE] = ode45(@cowell, [(ObserverOrbit.t(1)-gauss_ext(i).t2) 0], [gauss_ext(i).r2, gauss_ext(i).v2], options);
    % gauss_ext(i).r0 = statenew_GE(1,1:3);
    gauss_ext(i).coes_t0 = curtis_coe(gauss_ext(i).r0(1,:), gauss_ext(i).v0(1,:), mu);

    MonteCarlo_mod.r0_GE(i,:) = gauss_ext(i).r0;
    MonteCarlo_mod.v0_GE(i,:) = gauss_ext(i).v0;
    MonteCarlo_mod.a_t0_GE(i,:) = gauss_ext(i).coes_t0(end);
    MonteCarlo_mod.ecc_t0_GE(i,:) = gauss_ext(i).coes_t0(2);
    MonteCarlo_mod.inc_t0_GE(i,:) = gauss_ext(i).coes_t0(4);
    MonteCarlo_mod.RAAN_t0_GE(i,:) = gauss_ext(i).coes_t0(3);
    MonteCarlo_mod.AoP_t0_GE(i,:) = gauss_ext(i).coes_t0(5);
    MonteCarlo_mod.TA_t0_GE(i,:) = gauss_ext(i).coes_t0(6);


    % Gauss Extended ACO
    [gauss_ext_ACO(i).r0, gauss_ext_ACO(i).v0] = universal_variable(gauss_ext(i).r2, gauss_ext_ACO(i).v2, (ObserverOrbit.t(1)-gauss_ext(i).t2), mu, TOL);
    gauss_ext_ACO(i).coes_t0 = curtis_coe(gauss_ext_ACO(i).r0(1,:), gauss_ext_ACO(i).v0(1,:), mu);

    MonteCarlo_mod.r0_GE_ACO(i,:) = gauss_ext_ACO(i).r0;
    MonteCarlo_mod.v0_GE_ACO(i,:) = gauss_ext_ACO(i).v0;
    MonteCarlo_mod.a_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(end);
    MonteCarlo_mod.ecc_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(2);
    MonteCarlo_mod.inc_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(4);
    MonteCarlo_mod.RAAN_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(3);
    MonteCarlo_mod.AoP_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(5);
    MonteCarlo_mod.TA_t0_GE_ACO(i,:) = gauss_ext_ACO(i).coes_t0(6);


end

% Other Desirable info
MonteCarlo_mod.run_time_total = t_MC;
MonteCarlo_mod.arc_length_secs = 86400*(time_span(end)-time_span(1));
MonteCarlo_mod.error_arcseconds = noise; % noise added to data
MonteCarlo_mod.r_site_ECI_obs = Observation.r_site;
MonteCarlo_mod.JD_obs = time_span;
MonteCarlo_mod.obsAngles0 = Observation0.angles; % angles before MC sim (i.e. original first/last measurements without EXTRA noise added)
MonteCarlo_mod.r0_truth = TargetOrbit.r0;
MonteCarlo_mod.v0_truth = TargetOrbit.v0;
MonteCarlo_mod.obs_points = height(MonteCarlo_mod.r_site_ECI_obs);
if Lambert_Type == 1
    MonteCarlo_mod.lambert_solver = "Gauss";
elseif Lambert_Type == 2
    MonteCarlo_mod.lambert_solver = "UV (curtis)";
elseif Lambert_Type == 3
    MonteCarlo_mod.lambert_solver = "Gooding";
elseif Lambert_Type == 4
    MonteCarlo_mod.lambert_solver = "Izzo-Gooding";
elseif Lambert_Type == 5
    MonteCarlo_mod.lambert_solver = "Battin";
end


save(full, "MonteCarlo_mod")

% Display summary table
mean_or_median = 2; % display median
MontePlot_GA(MonteCarlo_mod, mean_or_median)





clearvars -except arc_length_counter noise_vect_counter arc_length_vect noise_vect