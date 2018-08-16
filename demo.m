% This is a script performming shadow removal algorithm
% described in the paper:
% A New Shadow Removal Method using Color-Lines. CAIP2017.
% If you use this code, please cite our paper.
%
% Author: Xiaoming Yu, 2017.

img_shadow=im2double(imread('./image/test.tif'));
[img_deshadow]=shadowRemoval(img_shadow,[]);

