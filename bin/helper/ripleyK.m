function K = ripleyK(locs,dist,box,method);
% RipleysK: Calculate K statistic
% 
% K = RipleysK(locs,dist, box,method) calculates G, the K statistic at each 
% distance for the data with x-y coordinates given by locs, and the
% bounding rectangle given by box=[minX maxX minY maxY].
% If method=0, no edge correction is applied.
% If method=1, points are only used if they are at least h units away from
% the edge.
%
% Note: The L statistic may be calculated from the K statistic as follows: 
%   L = sqrt(K/pi)-h;
%   

if nargin<4, method=1; end
[N,k] = size(locs);
if k~=2, error('locs must have two columns'); end
rbox = min([locs(:,1)'-box(1);box(2)-locs(:,1)';locs(:,2)'-box(3); box(4)-locs(:,2)'] );
% rbox is distance to box

DX = repmat(locs(:,1),1,N)-repmat(locs(:,1)',N,1);
DY = repmat(locs(:,2),1,N)-repmat(locs(:,2)',N,1);
DIST = sqrt(DX.^2+DY.^2);
DIST = sort(DIST);

if method==1
K = zeros(length(dist),1);
for k=1:length(K)
    I = find(rbox>dist(k));
    if length(I)>0
        K(k) = sum(sum(DIST(2:end,I)<dist(k)))/length(I);
    end
end
elseif method==0
    K = zeros(length(dist),1);
    for k=1:length(K)
        K(k) = sum(sum(DIST(2:end,:)<dist(k)))/N;
    end
end

lambda = N/((box(2)-box(1))*(box(4)-box(3)));
K = K/lambda;
