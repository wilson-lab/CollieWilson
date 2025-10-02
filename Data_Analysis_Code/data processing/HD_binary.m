function [acc_binary,exp_binary] = HD_binary(acc_panelpos,exp_panelpos,heading_range)
% fetch times when heading was within specified range
acc_binary = abs(acc_panelpos)<heading_range;
exp_binary = abs(exp_panelpos)<heading_range;
end