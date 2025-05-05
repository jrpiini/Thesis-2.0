%% This function implements Ansalone and Curti's GA using a Lambert solver with two optimization variables (rho_0, rho_f)
function [chromosomes, Elite] = GA_IOD_lambert_GEO(Observation, Lambert_Type, use_3_pts, ACO)
rng("shuffle")
tic

if use_3_pts % Eliminate all but first, middle, and last point
    tdif = Observation.t(end) - Observation.t(1);
    t_mid = tdif/2;
    [~,t2_ind] = min(abs(Observation.t-t_mid));
    Observation.r_site = [Observation.r_site(1,:); Observation.r_site(t2_ind,:); Observation.r_site(end,:)];
    Observation.angles = [Observation.angles(1,:); Observation.angles(t2_ind,:); Observation.angles(end,:)];
    Observation.LOS_measurements = [Observation.LOS_measurements(1,:); Observation.LOS_measurements(t2_ind,:); Observation.LOS_measurements(end,:)];
    Observation.t = [Observation.t(1); Observation.t(t2_ind); Observation.t(end)];
end
    

if Lambert_Type == 1
    lambert_solver = @lambert_Gauss;
elseif Lambert_Type == 2
    lambert_solver = @lambert_curtis_UV;
elseif Lambert_Type == 3
    lambert_solver = @lambert_gooding;
elseif Lambert_Type == 4
    lambert_solver = @IzzoGooding;
elseif Lambert_Type == 5
    lambert_solver = @lambert_battin;
end

mu = 398600.44189;
TOL = 1e-8;

max_generations = 150;

r_o_p = norm(Observation.r_site(1,:)); % position at perigee
v_o_p = sqrt(mu/r_o_p); % speed at perigee
a_T = 35786+6378+1000; % default
a_min = 20000; % default
% a_T = 35786+6378+1000; % modified
% a_min = 35786+6378-1000; % modified
max_v_T_perigee = ((2*mu)/r_o_p) * (1 - (r_o_p/(2*a_T)));
D = (v_o_p + max_v_T_perigee) * (Observation.t(end)-Observation.t(1));

ecc_max = (1 - (r_o_p/a_T));% ecc_max = .1;
% ecc_max = .05;

st_dev_threshold = 5; % 5 km, stops if st dev of ranges is less than this


%% New initialization

% Grid Preferences
num_points = 100000; % default value
rho0_min = 34000; % default value
rho0_max = 40000; % default value, max observable slant range for LEO (alt <= 2000), see max_slant_range.m

% Create initialization grid
range_points = initial_chromosomes(rho0_min, rho0_max, D, num_points);

% figure
% scatter(range_points(:,1),range_points(:,2))

% Apply restrictions to initial grid
count_chromosomes = 0;
for i = 1:height(range_points)
    rho_0 = range_points(i,1);
    rho_f = range_points(i,2);
    
    r0 = Observation.r_site(1,:) + (rho_0*Observation.LOS_measurements(1,:));
    rf = Observation.r_site(end,:) + (rho_f*Observation.LOS_measurements(end,:));
    [v0, vf] = lambert_solver(r0, rf, (Observation.t(end)-Observation.t(1)));

    if ACO
        v0 = calc_ACO(r0, v0); 
        vf = calc_ACO(rf, vf);
    end

    coe = curtis_coe(r0, v0, mu);
    a = coe(7);
    ecc = coe(2);

    if a <= a_T && a >= a_min && ecc >=0 && ecc <= ecc_max
        count_chromosomes = count_chromosomes + 1;
        chromosomes.ranges(count_chromosomes, 1:2) = [rho_0 rho_f];
        chromosomes.r0(count_chromosomes, 1:3) = r0;
        chromosomes.rf(count_chromosomes, 1:3) = rf;
        chromosomes.v0(count_chromosomes,:) = v0;
        chromosomes.vf(count_chromosomes,:) = vf;
    end

end

if count_chromosomes == 0
    error('Zero chromosomes produced from initialization.')
end

if count_chromosomes < 5
    warning('Less than 5 chromosomes produced from initialization. Recommend to increase grid points.')
end

if count_chromosomes > 1000 % Reduce population if too many candidates
    indices_kept = randperm(size(chromosomes.ranges, 1), 1000);
    chromosomes.ranges = chromosomes.ranges(indices_kept,:);
    chromosomes.r0 = chromosomes.r0(indices_kept,:);
    chromosomes.rf = chromosomes.rf(indices_kept,:);
    chromosomes.v0 = chromosomes.v0(indices_kept,:);
    chromosomes.vf = chromosomes.vf(indices_kept,:);
end

if count_chromosomes == 1
    st_dev = [10 10 10];
else
    st_dev = std(chromosomes.ranges);
end
while count_chromosomes < 1000
    % Apply mutation to get 1000 exactly
    selected = randi([1, count_chromosomes]);
    chromosomes.ranges(count_chromosomes+1,1) = chromosomes.ranges(selected, 1) + (st_dev(1)*randn);
    chromosomes.ranges(count_chromosomes+1,2) = chromosomes.ranges(selected, 2) + (st_dev(2)*randn);
    r0 = Observation.r_site(1,:) + (chromosomes.ranges(count_chromosomes+1,1)*Observation.LOS_measurements(1,:));
    rf = Observation.r_site(end,:) + (chromosomes.ranges(count_chromosomes+1,2)*Observation.LOS_measurements(end,:));
    [v0, vf] = lambert_solver(r0, rf, (Observation.t(end)-Observation.t(1)));
    if ACO
        v0 = calc_ACO(r0, v0); 
        vf = calc_ACO(rf, vf);
    end

    coe = curtis_coe(r0, v0, mu);
    a = coe(7);
    ecc = coe(2);

    if a <= a_T && a >= a_min && ecc >=0 && ecc <= ecc_max && chromosomes.ranges(count_chromosomes+1,1) >= rho0_min && chromosomes.ranges(count_chromosomes+1,2) >= rho0_min
        chromosomes.v0(count_chromosomes+1,:) = v0;   chromosomes.vf(count_chromosomes+1,:) = vf;
        chromosomes.r0(count_chromosomes+1, 1:3) = r0;
        chromosomes.rf(count_chromosomes+1, 1:3) = rf;
        count_chromosomes=count_chromosomes+1;
    end
end

% Store initial grid for regeneration use later
initials = chromosomes;
% chromosomes.ranges(1,:) = [10000 10];

% figure
% scatter(chromosomes.ranges(:,1), chromosomes.ranges(:,2))
% xlabel('\rho_0')
% ylabel('\rho_f')

clear range_points count_chromosomes

%% Begin iterating through generations
for i = 1:max_generations % begin GA

    % Fitness Evaluation
    for j = 1:height(chromosomes.ranges)
        
        prod_prev = 1;
        for k = 2:height(Observation.t)-1
            [r, ~] = universal_variable(chromosomes.r0(j,:), chromosomes.v0(j,:), Observation.t(k), mu, TOL); % Propagate chromosome orbit to obs times
            angles = rho2RaDec_topo(r - Observation.r_site(k,:)); % predicted RA, Dec for chromosome
            % LOS_predicted = (LOS_from_RADec(angles(1), angles(2))); % predicted LOS vector
            chromosomes.fitness(j) = dot(Observation.LOS_measurements(k,:), (LOS_from_RADec(angles(1), angles(2))))*prod_prev; % dot product of observation and prediction. Ideal value is 1
            % chromosomes.fitness(j) = dot(Observation.LOS_measurements(k,:), LOS_predicted) + prod_prev; 
            prod_prev = chromosomes.fitness(j);
        end

    end

    % Selection
    [~, sortedIndices] = sort(chromosomes.fitness, 'descend');

    % Choose 100 most elite for new generation
    for j = 1:100
        NewGen.ranges(j,:) = chromosomes.ranges(sortedIndices(j),:);
        NewGen.r0(j,:) = chromosomes.r0(sortedIndices(j),:);
        NewGen.rf(j,:) = chromosomes.rf(sortedIndices(j),:);
        NewGen.v0(j,:) = chromosomes.v0(sortedIndices(j),:);
        NewGen.vf(j,:) = chromosomes.vf(sortedIndices(j),:);        
    end

    clear chromosomes
    chromosomes.ranges = NewGen.ranges;
    chromosomes.r0 = NewGen.r0;
    chromosomes.rf = NewGen.rf;
    chromosomes.v0 = NewGen.v0;
    chromosomes.vf = NewGen.vf;
    clear NewGen

    % Mutation (to make up next 800 chromosomes for new generation)

    j = 1;
    st_dev = std(chromosomes.ranges);
    chromosomes.st_dev_range(i,:) = st_dev;
    while j <= 800
        selected = randi([1, 100]);
        % chromosomes.ranges(j+100,1) = chromosomes.ranges(selected, 1) + (10*randn);
        % chromosomes.ranges(j+100,2) = chromosomes.ranges(selected, 2) + (10*randn);
        % New
        chromosomes.ranges(j+100,1) = chromosomes.ranges(selected, 1) + (st_dev(1)*randn);
        chromosomes.ranges(j+100,2) = chromosomes.ranges(selected, 2) + (st_dev(2)*randn);

        r0 = Observation.r_site(1,:) + (chromosomes.ranges(j+100,1)*Observation.LOS_measurements(1,:));
        rf = Observation.r_site(end,:) + (chromosomes.ranges(j+100,2)*Observation.LOS_measurements(end,:));
        [v0, vf] = lambert_solver(r0, rf, (Observation.t(end)-Observation.t(1)));
        if ACO
            v0 = calc_ACO(r0, v0); 
            vf = calc_ACO(rf, vf);
        end

        
        coe = curtis_coe(r0, v0, mu);
        a = coe(7);
        ecc = coe(2);

        if ecc <= ecc_max && a <= a_T && a >= a_min
        chromosomes.v0(j+100,:) = v0;   chromosomes.vf(j+100,:) = vf;
        chromosomes.r0(j+100, 1:3) = r0;
        chromosomes.rf(j+100, 1:3) = rf;
        j=j+1;
        end

    end


    % Crossover (to make up next 50 chromosomes for new generation)
    counter = 1;
    while counter <= 50
        selected1 = randi([1, 100]);
        selected2 = randi([1, 100]);
        chromosomes.ranges(counter+900, 1) = chromosomes.ranges(selected2, 1);
        chromosomes.ranges(counter+900, 2) = chromosomes.ranges(selected1, 2);

        r0 = Observation.r_site(1,:) + (chromosomes.ranges(counter+900,1)*Observation.LOS_measurements(1,:));
        rf = Observation.r_site(end,:) + (chromosomes.ranges(counter+900,2)*Observation.LOS_measurements(end,:));
        [v0, vf] = lambert_solver(r0, rf, (Observation.t(end)-Observation.t(1)));
        if ACO
            v0 = calc_ACO(r0, v0); 
            vf = calc_ACO(rf, vf);
        end

        coe = curtis_coe(r0, v0, mu);
        a = coe(7);
        ecc = coe(2);

        if abs(chromosomes.ranges(counter+900, 1) - chromosomes.ranges(counter+900, 2)) < D && ecc <= ecc_max && a <= a_T && a >= a_min
            chromosomes.v0(counter+900,:) = v0;   chromosomes.vf(counter+900,:) = vf;
            chromosomes.r0(counter+900, 1:3) = r0;
            chromosomes.rf(counter+900, 1:3) = rf;
            counter = counter + 1;
        end

    end

    % Remaining 50 chromosomes are newly generated with same process as
    % initialization
    indices_newgen = randperm(size(initials.ranges, 1), 50);
    chromosomes.ranges(951:1000, :) = initials.ranges(indices_newgen, :);
    chromosomes.r0(951:1000, :) = initials.r0(indices_newgen, :);
    chromosomes.rf(951:1000, :) = initials.rf(indices_newgen, :);
    chromosomes.v0(951:1000, :) = initials.v0(indices_newgen, :);
    chromosomes.vf(951:1000, :) = initials.vf(indices_newgen, :); %% CHECK OTHER FUNCTIONS FOR THIS LATER!!!

    if st_dev(1) <= st_dev_threshold && st_dev(2) <= st_dev_threshold
        break
    end


end % end GA



% Final Sort using Fitness Function
% Fitness Evaluation
    for j = 1:height(chromosomes.ranges)

        prod_prev = 1;
        for k = 2:height(Observation.t)-1
            [r, v] = universal_variable(chromosomes.r0(j,:), chromosomes.v0(j,:), Observation.t(k), mu, TOL); % Propagate chromosome orbit to obs times
            angles = rho2RaDec_topo(r - Observation.r_site(k,:)); % predicted RA, Dec for chromosome
            LOS_predicted = (LOS_from_RADec(angles(1), angles(2))); % predicted LOS vector
            chromosomes.fitness(j) = dot(Observation.LOS_measurements(k,:), LOS_predicted)*prod_prev; % dot product of observation and prediction. Ideal value is 1
            % chromosomes.fitness(j) = dot(Observation.LOS_measurements(k,:), LOS_predicted) + prod_prev;
            prod_prev = chromosomes.fitness(j);
        end

    end

    % Selection of 100 elite
    [sortedFitness, sortedIndices] = sort(chromosomes.fitness, 'descend');
    for j = 1:100
        Elite.ranges(j,:) = chromosomes.ranges(sortedIndices(j),:);
        Elite.r0(j,:) = chromosomes.r0(sortedIndices(j),:);
        Elite.rf(j,:) = chromosomes.rf(sortedIndices(j),:);
        Elite.v0(j,:) = chromosomes.v0(sortedIndices(j),:);
        Elite.vf(j,:) = chromosomes.vf(sortedIndices(j),:);
        Elite.fitness(j) = chromosomes.fitness(sortedIndices(j));
    end

    Elite.st_dev = std(Elite.ranges);
% Elite.r0 = 0;
chromosomes.run_time = toc;


end