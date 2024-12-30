# ScarpOffsetAnalysis
Matlab code to calculate the vertical separation, dip slip, heave, and throw from a topographic profile across a fault scarp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      Christine Regalla %
%                                                      last modified     %
%                                                      3/23/16           %
%                                                                        %     
%                                                                        %
% This script will calcualte vertical separation, heave, throw, and      %
% fault slip from a topographic profile across a fault scarp. This       %
% code reads a tab demited text file of topographic profile data, which  %
% must be formatted in two columns of x (distance) and z (elevation)     %
% data points. Once text profile is entered, the profile can be saved as % 
% a .mat file and reloaded in subsequent runs. Midpoint and regressions  %
% through lower and upper surfaces can be chosen graphically or entered  %
% numerically. Requires Matlab Statistics Toolbox.                       %
%                                                                        %
% Inputs:                                                                %
%   filename.txt = tab delimited text file (without formatting)          %
%           containing x and z data for the topographic scarp profile    %
%                                                                        %
% Outputs:                                                               %
%    VS = vertical separation, calculated as the difference between the  %
%          elevation of the projection of upper and lower surfaces at    %
%           the scarp midpoint                                           %
                                                    %
%    ru = slope and intercept of regression through lower surface        %
%    rl = slope and intercept of regression through upper surface        %
%    ru_l = intersection point between lower fan surface and scarp face  %
%    rl_i = intersection point between upper fan surface and scarp face  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
