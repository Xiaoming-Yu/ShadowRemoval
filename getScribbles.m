function [ scribbles ] = getScribbles( img )
%Get scribbles to guide the shadow detection.
%
% The details of the algorithm are described in our paper: 
% A New Shadow Removal Method using Color-Lines. CAIP2017.
% If you use this code, please cite the paper.
%
%   Input arguments:
%   ----------------
%	img - A shadow image in the range [0,1], type: double
%
%   Output arguments:
%   ----------------
%   scribbles  - the scribbles selected by user
%
% Author: Xiaoming Yu, 2017.

   disp('Please select the lit area for shadow detection...');
   figure;imshow(img);
   [~,xs,ys] = freehanddraw(gca,'color','r','linewidth',3);
    close
    lightmap=zeros([size(img,1),size(img,2)]);
    for index=1:length(xs)
       lightmap(round(ys(index)),round(xs(index)))=1;
    end
   lightmap=logical(imdilate(lightmap,strel('sphere',3)));
   disp('Please select the shadow area for shadow detection');
   figure;imshow(img);
   [~,xs,ys] = freehanddraw(gca,'color','r','linewidth',3);
    close
    shadowmap=zeros([size(img,1),size(img,2)]);
    for index=1:length(xs)
       shadowmap(round(ys(index)),round(xs(index)))=1;
    end
   shadowmap=logical(imdilate(shadowmap,strel('sphere',3)));
   scribbles=lightmap+0.5*shadowmap;
end

