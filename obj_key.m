function f = obj_key(x)

    % --------- Load all required data (kept at the very top) ----------
    load Matrix&Target;               % -> Z1 (n x n x Nf), ff (1 x Nf, GHz)
    load num_feeds;         % -> num_feeds (scalar)

    % ------------------------ Basic params ----------------------------
    num_ports = num_of_ports;
    [~, ~, Nf] = size(S_origin); %#ok<ASGLU>
    Z0 = 50;
    freqGHz = ff(:).';              % 1 x Nf (GHz)

    % ---------------- Component value tables --------------------------
    Nc = numel(capacitor_valuemat);
    Ni = numel(inductor_valuemat);

    % ================= 1) Split decision vector x =====================
    idx = 1;
    x1_length = ceil(num_ports/4);
    x1  = x(idx : idx + x1_length - 1);            idx = idx + x1_length;
    x2  = x(idx : idx + num_feeds  - 1);           idx = idx + num_feeds;

    rem_len   = numel(x) - idx;                    % no x5
    num_loads = floor(rem_len/2);
    x3  = x(idx : idx + num_loads - 1);            idx = idx + num_loads;
    x4  = x(idx : idx + num_loads - 1);

    % ================= 2) x1 -> open/short bitmap =====================
    x1_norm   = (x1 + 1)/2 * 15;
    x1_bin    = dec2bin(round(x1_norm), 4);
    bits      = x1_bin(:)'; bits = bits(1:num_ports);
    internal_state = (bits == '1');                % 1=short, 0=open

    % ============== 3) Map feed indices (unique; avoid dup) ===========
    tmp2 = (x2 + 1)/2 * (num_ports-1) + 1;
    idx2 = min(max(round(tmp2),1), num_ports);
    feed_index = zeros(1, num_feeds);
    used = false(1, num_ports);
    for k = 1:num_feeds
        cand = idx2(k);
        if ~used(cand)
            sel = cand;
        else
            d = 1;
            while true
                low  = cand - d; high = cand + d;
                if low>=1 && ~used(low),  sel=low;  break; end
                if high<=num_ports && ~used(high), sel=high; break; end
                d = d + 1;
            end
        end
        feed_index(k) = sel;
        used(sel) = true;
    end

    % ============== 4) Map load indices (avoid feeds/used) ============
    tmp3 = (x3 + 1)/2 * (num_ports-1) + 1;
    idx3 = min(max(round(tmp3),1), num_ports);
    load_index = zeros(1, num_loads);
    for k = 1:num_loads
        cand = idx3(k);
        if ~used(cand)
            sel = cand;
        else
            d = 1;
            while true
                low  = cand - d; high = cand + d;
                if low>=1 && ~used(low),  sel=low;  break; end
                if high<=num_ports && ~used(high), sel=high; break; end
                d = d + 1;
            end
        end
        load_index(k) = sel;
        used(sel) = true;
    end

    % =========== 5) x4 -> component type & library indices ============
    value_index = x4(:).';
    load_sign   = sign(value_index);               % +1=L, -1=C
    elem_idx    = zeros(1, num_loads);

    vpos = find(load_sign < 0);                    % capacitors
    if ~isempty(vpos)
        tmpV = (abs(value_index(vpos))+1)/2 * (Nc-1) + 1;
        elem_idx(vpos) = min(max(round(tmpV),1), Nc);
    end
    ipos = find(load_sign >= 0);                   % inductors
    if ~isempty(ipos)
        tmpL = (value_index(ipos)+1)/2 * (Ni-1) + 1;
        elem_idx(ipos) = min(max(round(tmpL),1), Ni);
    end
    cap_vals = capacitor_valuemat(elem_idx(vpos));
    ind_vals = inductor_valuemat(elem_idx(ipos));

    % Per-frequency load impedances
    w = freqGHz * 2*pi * 1e9;
    Z_load = zeros(num_loads, Nf);
    if ~isempty(vpos), Z_load(vpos,:) = -1 ./ (1i * (cap_vals(:) * w)); end
    if ~isempty(ipos), Z_load(ipos,:) =  1i * (ind_vals(:) * w);        end

    % ================= 6) Equivalent S over feed ports =================
    S_tot = S_origin;
    internal_port_index = setdiff(1:num_ports, feed_index);
    Gamma_int = compute_internal_reflection(S_tot, internal_port_index, load_index, Z_load, internal_state, Z0, feed_index);
    SF = compute_equivalent_S(S_tot, feed_index, internal_port_index, Gamma_int);

    % ================= 7) Objectives ==================================
    P = numel(feed_index);
    match_thr  = 0.4;
    couple_thr = 0.3;
    pexp       = 3;

    mk = true(size(ff));                            % full band
    self_pen = zeros(P,1);
    for ii = 1:P
        Skk = abs(squeeze(SF(ii,ii,mk)));
        v   = max(0, Skk - match_thr) / (1 - match_thr);
        self_pen(ii) = mean(v.^pexp);
    end

    pair_vals = [];
    for ii = 1:P-1
        for jj = ii+1:P
            Sij = abs(squeeze(SF(ii,jj,mk)));
            v   = max(0, Sij - couple_thr) / (1 - couple_thr);
            pair_vals(end+1) = mean(v(:).^pexp); %#ok<AGROW>
        end
    end
    pair_mean = 0; if ~isempty(pair_vals), pair_mean = mean(sqrt(pair_vals)); end

    f = [ mean(self_pen);
          pair_mean;
          0 ] * 1000;
end

%% ----------------------- Local helper functions ------------------------
function Gamma_int = compute_internal_reflection(S, int_idx, l_idx, Z_load, state, Z0, feeding_port_index)
    m        = size(S,1);
    port_num = 1:m;
    Gamma_int = zeros(m, size(Z_load,2));
    for i = 1:m
        port = port_num(i);
        idx_load = find(l_idx == port, 1);
        if ~isempty(idx_load)
            ZL = Z_load(idx_load, :);
            Gamma_int(i, :) = (ZL - Z0) ./ (ZL + Z0);
        else
            Gamma_int(i, :) = (state(i) == 1) * -1 + (state(i) == 0) * +1;
        end
    end
    Gamma_int(feeding_port_index, :) = [];
end

function SF = compute_equivalent_S(S, feed_idx, int_idx, G_int)
    p  = numel(feed_idx);
    Nf = size(S,3);
    SF = zeros(p, p, Nf);
    for fi = 1:Nf
        Sfi = S(:, :, fi);
        A   = Sfi(feed_idx, feed_idx);
        B   = Sfi(feed_idx, int_idx);
        C   = Sfi(int_idx, feed_idx);
        D   = Sfi(int_idx, int_idx);
        M   = diag(1./G_int(:, fi)) - D;
        SF(:, :, fi) = A + B * (M \ C);
    end
end
