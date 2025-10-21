function avg_obj_speed = compute_avg_obj_speed(timebase, visobj_history)
% CREATED: 10/20/2025 - MC
sim_dur = size(visobj_history,2);
visobj_section = visobj_history(:,1:round(sim_dur/3));

dt = mean(diff(timebase));
pos_smooth = smoothdata(visobj_section, 1, 'movmean', 25);  % smooth
speed = abs(diff(pos_smooth) ./ dt);
avg_obj_speed = mean(speed(:), 'omitnan');
end
