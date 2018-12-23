function [offsetMap, cor_scribbles] = getOffset(img, cor_scribbles)
% Estimate offset of color-line
%
% The details of the algorithm are described in our paper: 
% A New Shadow Removal Method using Color-Lines. CAIP2017.
% If you use this code, please cite the paper.
%
% Input argument:  img - A shadow image in the range [0,1], type: double
%                  cor_scribbles  - A user scribbles for offset correction,
%                  in the range [0,1], type: double
% Output argument: offsetMap - offset map of the color-lines.
%
% Author: Xiaoming Yu, 2017. 
    if ~exist('cor_scribbles','var') || isempty(cor_scribbles)
        disp('Please select the shadow gradient region with similar material for offset correction...');
        figure;imshow(img);
        [~,xs,ys] = freehanddraw(gca,'color','r','linewidth',3);
        close
        guidemap=zeros([size(img,1),size(img,2)]);
       for index=1:length(xs)
           guidemap(round(ys(index)),round(xs(index)))=1;
       end
        guidemap=logical(imdilate(guidemap,strel('sphere',3)));
    else
        guidemap = cor_scribbles > 0.4;
    end
    r = img(:,:,1); g = img(:,:,2); b = img(:,:,3);
    selectData = cat(2, r(guidemap), g(guidemap), b(guidemap));
    [~,V]=pca(selectData);
    v=V';
    p=mean(selectData);
    ori=p-v*dot(p,v)/norm(v);
    for channel=1:3
        img(:,:,channel)=img(:,:,channel)-ori(channel);
    end
    imgnorm=sqrt(img(:,:,1).^2+img(:,:,2).^2+img(:,:,3).^2);
    offsetMap=zeros(size(img));
    v=v./norm(v);
    for channel=1:3
        offsetMap(:,:,channel)=img(:,:,channel).*ori(channel)./(v(channel).*imgnorm+eps);
    end
    
    function [d v] = pca(A)
        m = mean(A,1)'  ;
        for c=1:3
            A(:,c) = A(:,c) - m(c) ;
        end
        A = A' * A ;
        A = A / size(A,1) ;
        [V D] = eig(A) ;
        d = diag(D) ;
        v = V(:,3) ;
    end
end



       