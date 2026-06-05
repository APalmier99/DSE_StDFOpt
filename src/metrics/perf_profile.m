function hl = perf_profile(H, gate, logplot)
% Performance profile following Dolan & More' (2002)
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
      T(p,s) = nfevs;
    end
  end
end
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];
if nargin < 3, logplot = 0; end
r = T./repmat(min(T,[],2),1,ns);
max_ratio = max(max(r));
r(isnan(r)) = 2*max_ratio;
r = sort(r);
hl = zeros(ns,1);
for s = 1:ns
    [xs,ys] = stairs(r(:,s),(1:np)/np);
    if xs(1)==1
        vv = find(xs==1,1,'last'); xs = xs(vv:end); ys = ys(vv:end);
    end
    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = [char(lines(sl)) colors(sc) markers(sm)];
    if logplot
        hl(s) = semilogx(xs,ys,option1);
    else
        hl(s) = plot(xs,ys,option1);
    end
    hold on;
end
ax = ancestor(hl(1), 'axes'); ax.XAxis.Exponent = 0; xtickformat('%.0f');
if logplot
  axis([1 1.1*max_ratio 0 1]); twop = floor(log2(1.1*max_ratio)); set(gca,'XTick',2.^(0:twop))
else
  axis([1 1.1*max_ratio 0 1]);
end
