% NOTE:
%   'chromo' stores the decision variables and the corresponding objective
%   function values of the FINAL generation (the file 'chromo.mat' is
%   overwritten at each iteration, so it always holds the last generation).
%   The stored decision vectors can later be translated into topology
%   information via the 'obj_key' function.

clc; clear; close all;
global N M D PopCon gen
N = 300;                % population size
M = 3;                 % number of objectives (fixed to 3 as requested)
gen = 1000;           % generations

% --- reference points & init population ---
[Z, N] = UniformPoint(N, M);          % reference points for NSGA-III
[~, P] = funfun();                    % initial population: [N x D]

% --- evaluate population (constraints optional) ---
try
    [F, PopCon] = CalObj(P);          % F: [N x 3], PopCon<=0 feasible
catch
    F = CalObj(P); PopCon = [];
end


if isempty(PopCon)
    zmin = min(F, [], 1);
else
    feas = all(PopCon <= 0, 2);
    if any(feas), zmin = min(F(feas, :), [], 1);
    else,         zmin = min(F, [], 1);
    end
end

% --- evolution loop ---
for it = 1:gen
    % constraint violation (0 if no constraints)
    if isempty(PopCon)
        v = zeros(size(P,1), 1);
    else
        v = sum(max(0, PopCon), 2);
    end

    % selection -> variation
    idx = TournamentSelection(2, N, v);
    Q   = GA(P(idx, :));              % offspring

    % evaluate offspring
    try
        [Fq, Cq] = CalObj(Q); %#ok<NASGU>
    catch
        Fq = CalObj(Q); Cq = []; %#ok<NASGU>
    end

    % update ideal point & environmental selection
    zmin = min([zmin; Fq], [], 1);
    P    = EnvironmentalSelection([P; Q], N, Z, zmin);

    % re-evaluate survivors for next loop
    try
        [F, PopCon] = CalObj(P);
    catch
        F = CalObj(P); PopCon = [];
    end

    % --- plot (M=3) ---
    figure(1); clf;
    plot3(F(:,1), F(:,2), F(:,3), 'o', 'MarkerSize', 5); grid on; axis tight
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
    title(sprintf('Generation %d', it));

    % --- archive (decision + objectives) ---
    chromo = [P, F]; %#ok<NASGU>
    save('chromo.mat', 'chromo');

end
