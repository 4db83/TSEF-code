function z = lag(x,n,v)
% PURPOSE: creates a matrix or vector of lagged values:
% -------------------------------------------------------
% NOTE: works with timetable objects, but NOT with table.
% -------------------------------------------------------
% USAGE: z = lag(x,n,v)
% where: x = input matrix or vector, (nobs x k)
%        n = order of lag
%        v = (optional) initial values (default=NaN)
% e.g.
%     z = lag(x) creates a matrix (or vector) of x, lagged 1 observations
%         with initial value being set to zero. so needs to be trimmed to
%         get rid of initial 0.
%     z = lag(x,n) creates a matrix (or vector) of x, lagged n observations
%     z = lag(x,n,v) creates a matrix (or vector) of x, lagged n observations,
%         with initial values taking a value v.
% ------------------------------------------------------
% RETURNS: z = matrix (or vector) of lags (nobs x k)
% ------------------------------------------------------
% NOTES: if n <= 0, z = [] is returned. While you may find this
%        preverse, it is sometimes useful.
%-------------------------------------------------------
% SEE ALSO: mlag() 
%-------------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

switch(nargin)

case 1
   n = 1; v = NaN;
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];

case 2
   v = NaN;
   if n < 1
   z = [];
   else
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];
   end

case 3
   if n < 1
   z = [];
   else
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr(x,0,n)];
   end

otherwise
error('lag: wrong # of input arguments');
end

  
