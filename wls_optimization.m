function out = wls_optimization(in, data_weight, guidance, lambda)
% Energy function used in "A New Shadow Removal Method using Color-Lines". CAIP2017
%
%  Input arguments:
%  ----------------
%  in             - Input image (2-D, double, N-by-M matrix).   
%  data_weight    - High values indicate it is accurate, small values
%                   indicate it's not.
%  guidance       - Source image for the affinity matrix. Same dimensions
%                   as the input image IN. Default: log(IN)
%  lambda         - Balances between the data term and the smoothness
%                   term. Increasing lambda will produce smoother images.
%                   Default value is 0.05 
%
% This function is based on the implementation of the WLS Filer by Farbman,
% Fattal, Lischinski and Szeliski, "Edge-Preserving Decompositions for 
% Multi-Scale Tone and Detail Manipulation", ACM Transactions on Graphics, 2008
% The original function can be downloaded from: 
% http://www.cs.huji.ac.il/~danix/epd/wlsFilter.m
% Function is written by Berman, D. and Treibitz, T. and Avidan S.,
% and used in  "Non-Local Image Dehazing". CVPR2016.

small_num = 0.00001;

if ~exist('lambda','var') || isempty(lambda), lambda = 0.05; end

[h,w,isColor] = size(guidance);
k = h*w;
if isColor>1
guidance = mean(guidance,3);
end
if(size(lambda,1)==1)
    lambda=ones(size(data_weight)).*lambda;
end

% Compute affinities between adjacent pixels based on gradients of guidance
dy = diff(guidance, 1, 1);
dy = -lambda(1:end-1,:)./(sum(abs(dy).^2,3) + small_num);
dy = padarray(dy, [1 0], 'post');
dy = dy(:);

dx = diff(guidance, 1, 2); 
dx =-lambda(:,1:end-1)./(sum(abs(dx).^2,3) + small_num);
dx =padarray(dx, [0 1], 'post');
dx = dx(:);


% Construct a five-point spatially inhomogeneous Laplacian matrix
B = [dx, dy];
d = [-h,-1];
tmp = spdiags(B,d,k,k);

ea = dx;
we = padarray(dx, h, 'pre'); we = we(1:end-h);
so = dy;
no = padarray(dy, 1, 'pre'); no = no(1:end-1);

D = -(ea+we+so+no);
Asmoothness = tmp + tmp' + spdiags(D, 0, k, k);

% Normalize data weight
data_weight = data_weight - min(data_weight(:)) ;
data_weight = 1.*data_weight./(max(data_weight(:))+small_num);

% Make sure we have a boundary condition for the top line:
% It will be the minimum of the transmission in each column
% With reliability 0.8
reliability_mask = data_weight(1,:) < 0.6; % find missing boundary condition
in_row1 = min( in,[], 1);
data_weight(1,reliability_mask) = 0.8;
in(1,reliability_mask) = in_row1(reliability_mask);

Adata = spdiags(data_weight(:), 0, k, k);

A = Adata + Asmoothness;%单位整需要乘以一个权重再加上去，也就是说Adata=I.*w,所以前面构造拉普拉斯矩阵D时才没有像原始版本一样直接加上1，而是放在这里
b = Adata*in(:);%

% Solve
out = A\b;
% out = pcg(A,b,1e-4, 100, [], []);
out = reshape(out, h, w);
