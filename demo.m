% This is a script performming shadow removal algorithm
% described in the paper:
% A New Shadow Removal Method using Color-Lines. CAIP2017.
% If you use this code, please cite our paper.
%
% Author: Xiaoming Yu, 2017.

clear all;
img_shadow=im2double(imread('./image/test.tif'));
det_scribbles=im2double(imread('./image/det_scribbles.bmp'));
cor_scribbles=im2double(imread('./image/cor_scribbles.bmp'));
[img_deshadow]=shadowRemoval(img_shadow,det_scribbles,cor_scribbles);

