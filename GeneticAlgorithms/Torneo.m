function k = Torneo(aptitude)
    
    N = numel(aptitude);
    tao = round(N*0.3);

    Index = randi(N,[1 tao]);
    [~,i] = max(aptitude(Index));

    k = Index(i);