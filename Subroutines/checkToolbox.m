function [optim,symmath]=checkToolbox

v = ver;
optim = any(strcmp(cellstr(char(v.Name)), 'Optimization Toolbox'))
if optim==0
fprintf(1,'Optimization Toolbox is missing. Please install it first')
end
symmath = any(strcmp(cellstr(char(v.Name)), 'Symbolic Math Toolbox'))
if symmath==0
fprintf(1,'Symbolic Math Toolbox is missing. Please install it first')
end
