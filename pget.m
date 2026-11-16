function v = pget(P, name, default)
% PGET  Field of the parameter structure P, or a default when it is absent.
if isfield(P, name) && ~isempty(P.(name)), v = P.(name); else, v = default; end
end
