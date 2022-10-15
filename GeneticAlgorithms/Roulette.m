function k = Roulette (aptitude)

    aptitude_total = sum(aptitude); % get total of all aptitudes

    % p is a vector that contains the probablilty of each individual
    p = aptitude/aptitude_total;    
    
    r = rand();
    
    p_sum = 0;
    N = numel(aptitude);
    k = N;
    
    for i=1:N
        p_sum = p_sum + p(i);
    
        if p_sum >= r
            k = i;
            return
        end
    end