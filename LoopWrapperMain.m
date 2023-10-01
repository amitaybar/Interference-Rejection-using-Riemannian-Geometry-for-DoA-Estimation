close all
clc
clear
%%
% set(0,'DefaultFigureWindowStyle','docked');
%%

KFlags4Plot = 6;
for CurFlagInd=1:KFlags4Plot
    FlagsArray = zeros(1,KFlags4Plot);
    FlagsArray(CurFlagInd) = 1;
    
    LoopWrapper;
end

