function [right_ascension_mod,declination_mod,time_span_mod] = make_less_dense(right_ascension,declination,time_span,spacing)

% Makes data less dense according to spacing
% Ex: spacing = 1 means roughly 1 hz frequency measurements

JD_prev = time_span(1);
right_ascension_mod = right_ascension(1);
declination_mod = declination(1);
time_span_mod = time_span(1);

counter = 1;
for i = 1:height(time_span)
    t_dif = (time_span(i)-JD_prev)*86400;

    if t_dif >= spacing
        counter = counter + 1;
        right_ascension_mod(counter,:) = right_ascension(i);
        declination_mod(counter,:) = declination(i);
        time_span_mod(counter,:) = time_span(i);

        JD_prev = time_span(i);

    elseif i == height(time_span)
        counter = counter + 1;
        right_ascension_mod(counter,:) = right_ascension(i);
        declination_mod(counter,:) = declination(i);
        time_span_mod(counter,:) = time_span(i);
    end



end

