% print_specs.m  - print the computer specification for Section 3 of the paper
fprintf('MATLAB  : %s\n', version);
fprintf('OS      : %s\n', computer);
if ispc
    [~, cpu] = system('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor).Name"');  % Windows
    [~, mem] = memory; fprintf('RAM     : %.1f GB\n', mem.PhysicalMemory.Total/2^30);
elseif ismac
    [~, cpu] = system('sysctl -n machdep.cpu.brand_string');
    [~, ram] = system('sysctl -n hw.memsize'); fprintf('RAM     : %.1f GB\n', str2double(ram)/2^30);
else
    [~, cpu] = system('grep -m1 "model name" /proc/cpuinfo');
    [~, ram] = system('grep MemTotal /proc/meminfo'); fprintf('RAM     : %s', ram);
end
fprintf('CPU     : %s\n', strtrim(cpu));
