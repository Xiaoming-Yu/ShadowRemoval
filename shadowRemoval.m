function [img_deshaow] = shadowRemoval(img_shadow, scribbles)
% Perform shadow removal using color-lines.
%
% The details of the algorithm are described in our paper: 
% A New Shadow Removal Method using Color-Lines. CAIP2017.
% If you use this code, please cite the paper.
%
% Input arguments:
% img_shadow - A shadow image in the range [0,1], type: double
% scribbles  - A user scribbles gray image for shadow detection,
%              in the range [0,1], type: double
%
% Output argument:
% img_deshaow - shadow free image in the range [0,1], type: double.
%
% Author: Xiaoming Yu, 2017. 

% Validate input
[H,W,N_COLORS] = size(img_shadow);
if (N_COLORS ~= 3) % input verification
    error(['An RGB input image is required, while input ',...
        'has only ',num2str(N_COLORS),' dimensions']);
end
if ~exist('scribbles','var') || isempty(scribbles)
    [scribbles]=getScribbles(img_shadow);
else
end

% Offset correction
[offset]=getOffset(img_shadow);
img_c = zeros([H,W,3]);
if size(offset,1)==H
     img_c= img_shadow - offset;
else
    for color_idx=1:3
    img_c(:,:,color_idx) = img_shadow(:,:,color_idx) - offset(color_idx);
    end
end
% Finde color-lines and estimate luma
load('TR1000.mat');
Npoints = 1000;
ANG=TR.Points;
mdl = KDTreeSearcher(ANG);
% Original image clustering
img_1d=reshape(img_shadow,H*W,3);
radius = sqrt( img_shadow(:,:,1).^2 + img_shadow(:,:,2).^2 +img_shadow(:,:,3).^2 );
dist_norm = reshape(radius,[H*W,1]);
dist_unit_radius = bsxfun(@rdivide, img_1d, dist_norm);
color_ind = knnsearch(mdl, dist_unit_radius);
clu_ori = reshape(ANG(color_ind,:),H,W,3);
% Correction image clustering
img_c1d=reshape(img_c,[H*W,3]);
radius = sqrt( img_c(:,:,1).^2 + img_c(:,:,2).^2 +img_c(:,:,3).^2 );
dist_norm = reshape(radius,[H*W,1]);
dist_unit_radius = bsxfun(@rdivide, img_c1d, dist_norm);
color_ind = knnsearch(mdl, dist_unit_radius);
colorMax = accumarray(color_ind,radius(:),[Npoints,1],@max);
radius_max = reshape( colorMax(color_ind), H, W);
luma_estimation = (radius)./(radius_max+eps);
bin_count       = accumarray(color_ind,1,[Npoints,1]);
bin_count_map   = reshape(bin_count(color_ind),H,W);
bin_eval_fun    = @(x) min(1, x/50);
color_std = accumarray(color_ind,radius(:),[Npoints,1],@std);
radius_std = reshape( color_std(color_ind), H, W);
radius_eval_fun = @(r) min(1, 3*max(0.001, r-0.1));
radius_reliability = radius_eval_fun(radius_std./max(radius_std(:)));
data_term_weight   = bin_eval_fun(bin_count_map).*radius_reliability;
lambda = 0.001;
luma = wls_optimization(luma_estimation, data_term_weight, img_c, lambda);

% Shadow detection 
lightmap = scribbles>0.9; shadowmap = scribbles<0.6&scribbles>0.4;
limg=log(max(img_1d,eps));
trdata=[limg(lightmap(:),:);limg(shadowmap(:),:)];
trlabel=zeros([size(trdata,1),1]);
trlabel(1:sum(lightmap(:)))=1;
model = ClassificationKNN.fit(trdata,trlabel,'NumNeighbors',3);
shadow_ind = predict(model,limg);
shadow_ind=reshape(shadow_ind,H,W);
gsH = fspecial('gaussian',14,round(7));
shadow_ind = imfilter(shadow_ind,gsH,'replicate');
rmask=zeros([H,W]);
rmask(shadow_ind<0.5)=1;
rmask=bwareaopen(~bwareaopen(~rmask,100),200);
bd=getbd(rmask,ceil(max([H W])/512));
bl=bdspln(luma,rmask,bd,8);
bdv=bl.e-bl.s;
vw = arrayfun(@(x) norm(bdv(:,x)), 1:length(bl.s));
bd.w = zeros(1,length(bd.t)); bd.w(~bd.t) = vw;
pmask=getpmsk(bd,[H,W]);
umask = ~pmask&rmask;
lmask=~(umask|pmask);

% Umbra shadow removal
umask3d=reshape([umask umask umask],[H,W,3]);
lmask3d=reshape([lmask lmask lmask],[H,W,3]);
uind=color_ind;
lind=color_ind;
uind(umask(:)==0)=Npoints+1;
lind(lmask(:)==0)=Npoints+1;
lstd=zeros([Npoints+1 3]);
ustd=zeros([Npoints+1 3]);
lmean=zeros([Npoints+1 3]);
umean=zeros([Npoints+1 3]);
for channel=1:3
lstd(:,channel)=accumarray(lind,img_1d(:,channel),[Npoints+1 1],@std);
ustd(:,channel)=accumarray(uind,img_1d(:,channel),[Npoints+1 1],@std);
lmean(:,channel)=accumarray(lind,img_1d(:,channel),[Npoints+1 1],@mean);
umean(:,channel)=accumarray(uind,img_1d(:,channel),[Npoints+1 1],@mean);
end
% Use the average statistis to relight the umbra region that has no reference lit pixels
sMaterialindex=sum(umean,2)>0;
NRMaterialindex=sMaterialindex & (sum(lmean,2)==0);
if find(NRMaterialindex)
    sNormalMaterialindex=sMaterialindex & ~NRMaterialindex;
    meanscale=mean(lmean(sNormalMaterialindex,:)./umean(sNormalMaterialindex,:),1);
    stdscale=mean(lstd(sNormalMaterialindex,:)./ustd(sNormalMaterialindex,:),1);
    for index=1:3
        lmean(NRMaterialindex,index)=umean(NRMaterialindex,index)* meanscale(index);
        lstd(NRMaterialindex,index)=ustd(NRMaterialindex,index)* stdscale(index);
    end
end
% Relight umbra region.
std_scale=lstd./max((ustd+eps),0.01);
udeimg=(img_1d-umean(uind,:)).*std_scale(uind,:)+lmean(uind,:);

% Penumbra shadow removal
pmask3d=reshape([pmask pmask pmask],[H,W,3]);
pind=color_ind;
pind(pmask(:)==0)=Npoints+1;
uluma_mean=accumarray(uind,luma(:),[Npoints+1 1],@mean);
lluma_mean=accumarray(lind,luma(:),[Npoints+1 1],@mean);
uscale=img_1d./udeimg;
uscale_mean=zeros([Npoints+1 3]);
pscale=zeros([H*W,3]);
for channel=1:3
    uscale_mean(:,channel)=accumarray(uind,uscale(:,channel),[Npoints+1 1],@mean);
    pscale(:,channel)=uscale_mean(pind,channel)+(1-uscale_mean(pind,channel)).*(luma(:)-uluma_mean(pind))./(lluma_mean(pind)-uluma_mean(pind));
end
pscale=reshape(pscale,H,W,3);
uscale=reshape(uscale,H,W,3);
pscale(pmask3d==0)=0;
uscale(umask3d==0)=0;
scale_estimatiom=pscale+uscale;
scale_estimatiom(lmask3d==1)=1;
img_deshaow=img_shadow./scale_estimatiom;

% Illumination Optimization
anomalousValue = (((img_deshaow(:,:,1)>=1 | img_deshaow(:,:,2)>=1 | img_deshaow(:,:,3)>=1)...
        |(img_deshaow(:,:,1)<0.01 | img_deshaow(:,:,2)<0.01 | img_deshaow(:,:,3)<0.1))&logical(umask));
scale_estimatiom(scale_estimatiom>1)=1;
scale_estimatiom(scale_estimatiom<0)=1;
scale_estimatiom(isnan(scale_estimatiom))=1;
scale = zeros([H,W,3]);
bp=boundarymask(pmask);
data_term_weight(bp ==1)=data_term_weight(bp==1).*0.4;
data_term_weight(anomalousValue)=0.4;
data_term_weight(lmask==1)=1;
lambda = lmask*0.001+pmask*0.001+umask*0.001;

for color_idx=1:3
    scale(:,:,color_idx)=wls_optimization(scale_estimatiom(:,:,color_idx),data_term_weight,luma,lambda);
end
img_deshaow=img_shadow./scale;

figure,
subplot(2,4,1);imshow(img_shadow);title('Input');
subplot(2,4,2);imshow(clu_ori);title('Original image clustering');
subplot(2,4,3);imshow(luma_estimation);title('Initial luma estimation');
subplot(2,4,4);imshow(lmask*0.3+pmask*0.7+umask);title('Shadow mask');

subplot(2,4,5);imshow(img_c);title('Correction image');
subplot(2,4,6);imshow(reshape(ANG(color_ind,:),H,W,3));title('Correction image clustering');
subplot(2,4,7);imshow(luma);title('Refined luma estimation');
subplot(2,4,8);imshow(img_deshaow);title('Shadow free image'); 
end

