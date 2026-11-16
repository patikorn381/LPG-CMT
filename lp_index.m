function s = lp_index(m)
% Radial index of an LP0m label: '2'..'9', and ',10', ',11', ... from 10 on
% (LP02 ... LP09, LP0,10 ...).
if m < 10, s = sprintf('%d', m); else, s = sprintf(',%d', m); end
end
