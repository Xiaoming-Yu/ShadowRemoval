function [prolen] = get2Dprolen(vec,vecpro)
% [prolen] = get2DproLen(vec,vecpro) This function gets the projection length
% 
% INPUTS:
%   vec: a 2D vector to be projected
%   vecpro: 2 2D vector to be projected on
% OUTPUTS:
%   prolen: length of projection
% 
% Function is written by Han Gong, University of Bath, UK 

prolen = dot(vec,vecpro)/norm(vecpro);