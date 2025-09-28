function PopObj = CalObj(PopDec)
    global M
    NN = size(PopDec,1);
    % ---- loads at the very top ----
    load num_feeds
    PopObj = zeros(NN, M);
    for k = 1:NN
        f = obj_key(PopDec(k,:));
        PopObj(k,:) = f(:).';   % ensure row vector
    end
end
