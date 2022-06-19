
function [D,PC] = kmeans_opt(X,ToTest,distfun)
%%% [D,PC]=kmeans_opt(X,MAX) returns the cluster distance and explained variance for each datapoint in
%%% vector X.

%%% sebastien.delandtsheer@uni.lu
%%% sebdelandtsheer@gmail.com
%%% Thomas.sauter@uni.lu
warning('off','all')

% unit-normalize
X = X(~isnan(X(:,1)),:);
MIN=min(X); MAX=max(X); 
X=(X-MIN)./(MAX-MIN);
%tot_withinss = [sum(d**2) for d in dist]  # Total within-cluster sum of squares
%totss = sum(pdist(X)**2)/X.shape[0]       # The total sum of squares
%betweenss = totss - tot_withinss          # The between-cluster sum of squares

[~,S] = pca(X);
tot_witSS = zeros(size(ToTest)); % initialize the WSS
tot_SS = sum(pdist(S','correlation').^2);
for c = 1:numel(ToTest) % for each sample
    [~,~,sumd,d] = kmeans(S,ToTest(c),'Distance',distfun); %compute the sum of intra-cluster distances
    tot_witSS(c) = sum(min(d,[],2).^2);
    c
end
bet_SS = (tot_SS-tot_witSS)/tot_SS;

Var = -diff(D); % Calculate distance improvement for each step
PC = cumsum(Var)/(D(1)-D(end)); % Calculate the contribution to overall variance for each step
end