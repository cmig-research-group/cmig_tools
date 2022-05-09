function state = mmil_check_nargs(nargs, expected_nargs, errorMsg)
%function mmil_check_nargs(nargs, expected_nargs, errorMsg)
%
% Purpose: to have a common function to check and handle 
%   this simple parameter passing check.
%
%   Currently, the function compares the # of args and, 
%   if they are different, throws a formatted error message.
%
%   NOTE that it is not important to be careful with logging
%   here, as this is either a programming error or one that happens
%   at the start of script execution.
%
% Output:
%  state: 0 for failure (nargs don't match), 1 for success (args match)
%
% Created By:       Ben Cipollini on 08/15/2007
% Last Modified By: Ben Cipollini on 08/20/2007

  if (nargs >= expected_nargs)
    state = 1;
    
  else
    state = 0;
      
    all_callers = dbstack;
    caller  = all_callers(2).name;
	
  	if (~exist('errorMsg', 'var'))
	    errorMsg = sprintf('There are %d required argument(s) to call %s', ...
	                       expected_nargs, caller);
    end;
	
    % dump the help
    fprintf('\n');
    help(caller);
    fprintf('\n');
    
    % dump the error
    if (nargs ~= 0)
      error('\t%s', errorMsg );
    end;
  end;
