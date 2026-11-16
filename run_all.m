% RUN_ALL  Regenerate every table and figure of the paper and write a log
% (results/summary.txt) that can be sent back for checking.
if ~exist('results','dir'), mkdir('results'); end
if exist(fullfile('results','summary.txt'),'file'), delete(fullfile('results','summary.txt')); end
diary(fullfile('results','summary.txt')); diary on;
setappdata(0, 'lpg_t0', tic);   % survives the 'clear' inside the scripts
fprintf('==== print_specs ====\n');        print_specs;
fprintf('\n==== lpg_main ====\n');         lpg_main;
fprintf('\n==== validation_specs ====\n'); validation_specs;
fprintf('\n==== compare_singh ====\n');    compare_singh;
fprintf('\n==== dispersion_check ====\n'); dispersion_check;
fprintf('\n==== convergence_check ====\n'); convergence_check;
fprintf('\n==== tolerance_analysis ====\n'); tolerance_analysis;
fprintf('\nALL DONE in %.1f s. Next, in this folder: python verify_independent.py   then   python vector_check.py\n', toc(getappdata(0, 'lpg_t0')));
diary off;
