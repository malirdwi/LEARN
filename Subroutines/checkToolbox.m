function checkToolbox

v = ver;
optim = any(strcmp(cellstr(char(v.Name)), 'Optimization Toolbox'));
if optim==0
error('Optimization Toolbox is missing. Please install it before running this subroutine');
elseif ~exist('linprog','file')
        error('Optimization Toolbox is missing. Please install it before running this subroutine');
    end
symmath = any(strcmp(cellstr(char(v.Name)), 'Symbolic Math Toolbox'));
if symmath==0
error('Symbolic Math Toolbox is missing. Please install it before running this subroutine');
elseif ~exist('syms','file')
        error('Symbolic Math Toolbox is missing. Please install it before running this subroutine');
    end
       
    
if ~exist('cvx_begin','file')
        error('cvx is missing. Please install it via cvxr.com/cvx/download before running this subroutine');
    end
