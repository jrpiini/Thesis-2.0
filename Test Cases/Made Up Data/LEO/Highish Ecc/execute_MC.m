clear; close all; clc

arc_length_vect = [15 30 60]; % seconds
noise_vect = [3 10]; % arcseconds

for arc_length_counter = 1:length(arc_length_vect)
    % if arc_length_vect(arc_length_counter) == 30
    %     noise_vect = 2;
    % else
    %     noise_vect = [1 2 3];
    % end
	for noise_vect_counter = 1:length(noise_vect)
		arc_length = arc_length_vect(arc_length_counter);
		noise = noise_vect(noise_vect_counter);
		run init_Lambert_GA_6_58_MOD.m
	end
end