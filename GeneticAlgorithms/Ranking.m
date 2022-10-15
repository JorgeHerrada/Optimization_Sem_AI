function k = Ranking(aptitude)

    [~,Index] = sort(aptitude,'descend');
    N = numel(aptitude);

    rank = N:-1:1;
    rank_total = sum(rank);

    r = rand();
    p_sum = 0;
    k = N;

    for i=1:N
        p_sum = p_sum + rank(i)/rank_total;

        if p_sum>=r
            k = Index(i);
            return
        end
    end

