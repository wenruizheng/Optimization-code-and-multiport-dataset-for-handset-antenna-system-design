function [PopObj, PopDec] = funfun()
    global M D lower upper encoding N PopCon cons cons_flag name;

    f_num = M;

    % ---------- Load required data at the very top ----------
    load Matrix&Target
%     load freindex_mat

    % ---------- Decision vector layout parameters ----------
    num_feeds = 4;
    num_loads = 14;

    % x = [x1 | x2 | x3 | x4]
    % x1 packs open/short bits (4 bits per port entry)
    x1_length = ceil(num_of_ports / 4);
    x_num = x1_length + num_feeds + 2 * num_loads;   % Decision Variation code

    % presets
    feed_pointed = [1 50];
    load_pointed = [51 52 53 54];

    % component libraries
    capacitor_valuemat = -[0.5:0.1:2 2.2 2.4 2.7 3 3.3 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1 10 12 15 18 22 27 33 39 47 68] * 1e-12;
    inductor_valuemat  =  [0.6:0.1:2.4 2.7 2.9 3.0 3.3 3.4 3.6 3.9 4.3 4.7 5.1 5.6 6.2 6.8 7.5 8.2 9.1 10 11 12 15 18 33] * 1e-9;

    % persist parameters for obj_key (single MAT as before)
    save('num_feeds.mat', 'feed_pointed', 'load_pointed', ...
         'inductor_valuemat', 'capacitor_valuemat', ...
         'f_num', 'x_num', 'num_feeds', 'num_loads');

    % ---------- Population initialization ----------
    D = x_num;
    lower = -ones(1, D);
    upper =  ones(1, D);
    encoding = 'real';

    switch encoding
        case 'binary'
            PopDec = randi([0, 1], N, D);
        case 'permutation'
            [~, PopDec] = sort(rand(N, D), 2);
        otherwise
            PopDec = unifrnd(repmat(lower, N, 1), repmat(upper, N, 1));
    end

    % ---------- Constraints (kept as-is) ----------
    cons = zeros(size(PopDec, 1), 1);
    cons_flag = 1;
    PopCon = cons;

    % ---------- Evaluate objectives ----------
    PopObj = zeros(N, M);
            for KKK = 1:N
                x = PopDec(KKK, :);
                f = obj_key(x);
                tempObj = zeros(1, length(f));
                for k_i = 1:length(f)
                    tempObj(k_i) = f(k_i);
                end
                PopObj(KKK, :) = tempObj;
            end

end
