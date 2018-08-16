function bd = getbd(msk,fsp)
%GETBD obtains the boundary of shadow
%
% Function is written by Han Gong, University of Bath, UK 

imhw = size(msk); % size of image
% get boundaries of sampling areas
[bdp,~,N] = bwboundaries(msk); % boundary
bdl = cellfun(@(x) size(x,1),bdp); % length of each boundary
pl = find(bdl>3); % remove too short boundaries
bdpf = bdp(pl); % filtered boundaries
bdnf = cell(size(bdpf)); % boundary normal
bdtf = cell(size(bdpf)); % type of boundary point
bdmf = cell(size(bdpf)); % maker of boudnary number

for i = 1:length(pl)
    cbd = flipud(bdpf{i}'); ptlen = size(cbd,2); % current bondary
    % compute boundary normal
    t = cumsum([0,sqrt((sum(diff(cbd,1,2).^2)))]); % x^2+y^2 计算边界长度
    cv = csaps(t,cbd,1/(1 + fsp^3/0.6)); %拟合出长度与位置的关系
    v = fntlr(cv,2,t); % derivatie vector 泰勒2阶系数，即是该点数据以及一阶导
    R = [0,1;-1,0]; if pl(i)>N, R = -R; end % rotation matrix
    cbdn = R*v(3:4,:); % normal vector
    cbdn = bsxfun(@rdivide,cbdn,sqrt(sum(cbdn.^2)));
    % exclude boudnary at image boarder
    bm = cbd(1,:)==1|cbd(1,:)==imhw(2)|cbd(2,:)==1|cbd(2,:)==imhw(1);
    cdt = zeros(ptlen,1); cdt(bm) = -1; % mark points at image boarder
    % curature sampling
    npt = fnval(cv,t); ksp = 20;
    curv = cumsum(sat(abs(LineCurvature2D(npt')),0,ksp)); curv = curv/ksp;
    clidx = [1;find(diff(fix(curv),1)>0);ptlen]; % selected bounday index
    lidx = union(transpose(1:fsp:ptlen),clidx); % combined boundary index
    bdpf{i} = transpose(cbd(:,lidx)); bdnf{i} = transpose(cbdn(:,lidx));
    bdtf{i} = cdt(lidx); bdmf{i} = ones(numel(lidx),1)*i;
end

bd.p = transpose(cell2mat(bdpf)); bd.n = transpose(cell2mat(bdnf));
bd.t = cell2mat(bdtf); bd.m = cell2mat(bdmf);

end
