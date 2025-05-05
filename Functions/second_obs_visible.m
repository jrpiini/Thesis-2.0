function [vis, RA_error, Dec_error] = second_obs_visible(r0_truth, v0_truth, r0_IOD, v0_IOD, T_truth, T_error, T_IOD, r_site, r_site_times, JD0)

if isnan(T_IOD)
    vis = 0;
    return
end

options = odeset('AbsTol',1e-8, 'RelTol',1e-8);

FOV = 0.7;

% % Testing
[~,lower_bound] = min(abs(r_site_times - (JD0 + T_truth/86400 - T_error/86400)));
[~,upper_bound] = min(abs(r_site_times - (JD0 + T_truth/86400 + T_error/86400)));
[~,mid] = min(abs(r_site_times - (JD0 + T_IOD/86400)));

r_site_IOD = r_site(mid,:);


time_vect = r_site_times(lower_bound:upper_bound);
r_site = r_site(lower_bound:upper_bound,:);

if isscalar(time_vect)
    vis = 0;
    return
end

[~, statenew_IOD] = ode45(@cowell, [0 (r_site_times(mid)-JD0)*86400], [r0_IOD v0_IOD], options);
r_IOD = statenew_IOD(end,1:3);
rho_IOD = r_IOD - r_site_IOD;


[~, stateout_truth] = ode45(@cowell, [0 86400*(time_vect(1)-JD0)], [r0_truth v0_truth], options);
[~, statenew_truth] = ode45(@cowell, 86400*(time_vect-time_vect(1)), [stateout_truth(end,:)], options);


vis = 0;
RA_error = zeros(height(time_vect),1);
Dec_error = RA_error;
for i = 1:length(time_vect)
    
    r_truth = statenew_truth(i,1:3);
    rho_truth = r_truth - r_site(i,:);
    angles = rho2RaDec_topo(rho_truth);
    RA_truth = angles(1); Dec_truth = angles(2);
    angles = rho2RaDec_topo(rho_IOD);
    RA_IOD = angles(1); Dec_IOD = angles(2);

    RA_error(i) = abs(RA_IOD - RA_truth);
    Dec_error(i) = abs(Dec_IOD - Dec_truth);

    if RA_error(i) <= (FOV/2) && Dec_error(i) <= (FOV/2)
        vis = 1;
        return
    end

end


end