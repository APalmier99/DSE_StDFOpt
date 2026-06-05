function hl = data_profile(H,N,gate)
% Data profile following More' & Wild (2009)
[nf,np,ns] = size(H);
for j = 1:ns
    for i = 2:nf
      H(i,:,j) = min(H(i,:,j),H(i-1,:,j));
    end
end
prob_min = min(min(H),[],3);
prob_max = H(1,:,1);
T = zeros(np,ns);
for p = 1:np
    cutoff = prob_min(p) + gate*(prob_max(p) - prob_min(p));
    for s = 1:ns
        nfevs = find(H(:,p,s) <= cutoff,1);
        if isempty(nfevs)
            T(p,s) = NaN;
        else
            T(p,s) = nfevs/N(p);
        end
    end
end
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
max_data = max(max(T));
T(isnan(T)) = 2*max_data;
T = sort(T);
hl = zeros(ns,1);
for s = 1:ns
    [xs,ys] = stairs(T(:,s),(1:np)/np);
    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];
    hl(s) = plot(xs,ys,option1);
    hold on;
end
axis([0 1.1*max_data 0 1]);
