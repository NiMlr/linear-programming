% Last modified: 12.07.2017

function [optimum, utility] = sensitivityAnalysis(A, b, c, base, modification1, modification2)
% This is a function that implements the sensitivity analysis.
%
% Note that any solutions where xi < 0 are disregarded and that multiple 
% changes to the problem can not be made in one step. 
%
% For example, a new variable that results in modification of A and c
% counts as a single problem change. When multiple modifications are
% necessary please set the modification1 variable to the modified A.

% Input:
%   A - nConditions x nVariables
%       The initial matrix (normal form) with the coeffitients 
%       (usually left side) of the boundary condition-equations
%       (along the lines). Note, that the equations must be of "<="-type.
%   b - nConditions x 1
%       The initial limits vector of the boundary conditions
%       (usually right side).
%   c - 1 x nVariables
%       The initial coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
%   modification1 - size of the first modified component
%       The first modified component of changed 
%       optimization problem (normal form).
%       For single component change this can be A,b or c.
%       For multiple component change this must be A.
%   modification2 (optional) - size of the second modified component
%       The second modified component of changed 
%       optimization problem (normal form).
%       Note that this input is only required if more than one of the
%       components (A,b,c) was modified. This input can not be
%       the modified A.
%
% Output:
%   optimum - nVariables x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
% ------------------------------------------------------------------------
% Example input:
%   Conditions: 6 x1 + 15x2 <= 4500        [ 6 15;
%               4 x1 + 5 x2 <= 2000 -> A =   4  5;
%               20x1 + 10x2 <= 8000          20 10 ]
%                               |             [ 4500;
%                                --> limits =   2000;
%                                               8000 ];
%   Utility-function: max 16x1 + 32x2  ->  utilityCoeffitients = [ 16 32 ];
%   The optimal starting base for this problem: startingBase = [ 1 2 5 ];
%   New Condition:
%       4 x1 + x2 <= 1500   -> modification1 = [A;[4 1]]
%                      |
%                       --> modification2 = [limits; 1500];
% ------------------------------------------------------------------------

[A, c] = normal2standard(A, c);                                              % switch to standard form

if nargin == 5                                                             % check how many modifications were made
    switch detectModificationType(modification1)
        case 1                                                             % case where only values (no dimensions) of A were modified
            [optimum,utility] = modifiedConditions(modification1,b,c,base);
        case 2                                                             % case where only values (no dimensions) of b were modified
            [optimum,utility] = modifiedLimits(A,modification1,c,base);
        case 3                                                             % case where only values (no dimensions) of c were modified
            [~,modification1] = normal2standard(A,modification1);
            [optimum,utility] = modifiedCost(A,b,modification1,base);
        otherwise
            inputError();return;
    end
else
    switch detectModificationType(modification1,modification2)
        case 1                                                             % case where a new variable was introduced (A and c were modifed)
            [modifiedA,modifiedc] = normal2standard(modification1,...
                modification2);
            base = correctBase(size(A,2),base);
            [optimum,utility] = newVariable(modifiedA,b,modifiedc,base);
        case 2                                                             % case where a new condition was introduced (A and b were modifed)
            modifiedA = normal2standard(modification1);
            [optimum,utility] = newCondition(modifiedA,modification2,c,...
                base);
        otherwise
            inputError();return;
    end

end
end



function indicator = detectModificationType(modification1, modification2)
% Function that detects the different types of modifications.
%
% Input:
%   modification1 - size of the first modified component
%       The first modified component of changed 
%       optimization problem (normal form).
%       For single component change this can be A,b or c.
%       For multiple component change this must be A.
%   modification2 (optional) - size of the second modified component
%       The second modified component of changed 
%       optimization problem (normal form).
%       Note that this input is only required if more than one of the
%       components (A,b,c) was modified. This input can not be
%       the modified A.
%
% Output:
%   indicator - scalar
%       A value that indicates which modification has been made.
sizeOfModification1 = size(modification1);
if nargin==1                                                               % single modification
    if length(sizeOfModification1)==2
        if sizeOfModification1(1)==1                                       % row vector => change in c
            indicator = 3;
        elseif sizeOfModification1(2)==1                                   % column vector => change in b
            indicator = 2;
        else
            indicator = 1;                                                 % change in A
        end
    else
        indicator = nan;
    end
else                                                                       % multiple modifications
    sizeOfModification2 = size(modification2);
    if sizeOfModification2(1)==1                                           % row vector => new variable
            indicator = 1;
    elseif sizeOfModification2(2)==1                                       % column vector => new condition
            indicator = 2;
    else
            indicator = nan;
    end
end
end

function [optimum, utility] = modifiedConditions(modifiedA, b, c, base)
% Function that implements the routine for modified condition 
% values in A.
%
% Input:
%   modifiedA - nConditions x nVariables
%       The matrix (normal form) with the modified coeffitients 
%       (usually left side) of the boundary condition-equations
%       (along the lines). Note, that the equations must be of "<="-type.
%   b - nConditions x 1
%       The initial limits vector of the boundary conditions
%       (usually right side).
%   c - 1 x nVariables
%       The initial coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
% Output:
%   optimum - nVariables x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
modifiedA = normal2standard(modifiedA);
complement = getComplement(base,size(modifiedA,2));
if invertibleAndValid(modifiedA(:,base),b)                                 % if valid solution
    complementCost = c(complement)-(c(base)*inv(modifiedA(:,base))...      
        *modifiedA(:,complement));
    if max(complementCost)>0                                               % check optimality
        baseSuboptimal();optimum = nan;utility = nan;
    else
        [optimum,utility] = computeOptimum(modifiedA, b, c, base);
    end
end
end

function [optimum, utility] = modifiedLimits(A, modifiedb, c, base)
% Function that implements the routine for modified limits in b.
%
% Input:
%   A - nConditions x nVariables
%       The initial matrix (standard form) with the coeffitients 
%       (usually left side) of the boundary condition-equations
%       (along the lines). Note, that the equations must be of "<="-type.
%   modifiedb - nConditions x 1
%       The modified limits vector of the boundary conditions
%       (usually right side). Note that only value (not dimension) 
%       modifications are handled.
%   c - 1 x nVariables
%       The initial coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
% Output:
%   optimum - nVariables x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
if invertibleAndValid(A(:,base),modifiedb)                                 % if still valid
    [optimum,utility] =  computeOptimum(A,modifiedb,c,base);
else
    solutionInvalid();optimum = nan;utility = nan;
end
end

function [optimum, utility] = modifiedCost(A, b, modifiedc, base)
% Function that implements the routine for modified cost in c.
%
% Input:
%   A - nConditions x nVariables
%       The initial matrix (standard form) with the coeffitients 
%       (usually left side) of the boundary condition-equations
%       (along the lines). Note, that the equations must be of "<="-type.
%   b - nConditions x 1
%       The initial limits vector of the boundary conditions
%       (usually right side).
%   modifiedc - 1 x nVariables
%       The modified coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm and that this function
%       only handles value (not dimension) modifications.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
% Output:
%   optimum - nVariables x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
complement = getComplement(base,size(A,2));
complementCost = modifiedc(complement)-(modifiedc(base)*inv(A(:,base))*... % "shadow" cost of the non-base variables
    A(:,complement));
if max(complementCost)>0                                                   % check optimality
    baseSuboptimal();optimum = nan;utility = nan;
else
    [optimum,utility] = computeOptimum(A, b, modifiedc, base);
end
end

function [optimum, utility] = newVariable(modifiedA, b, modifiedc, base)
% Function that implements the routine for an aditional variable.
%
% Input:
%   modifiedA - nConditions x nVariables+1
%       The matrix (normal form) with the appended coefficients for the new
%       variable. Note, that the equations must be of "<="-type.
%   b - nConditions x 1
%       The initial limits vector of the boundary conditions
%       (usually right side).
%   modifiedc - 1 x nVariables+1
%       The modified coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
% Output:
%   optimum - nVariables+1 x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
newVarIndex = length(modifiedc)-length(base);
rNew = modifiedc(newVarIndex)-(modifiedc(base)*inv(modifiedA(:,base))*...  % compute new variable "shadow" cost
    modifiedA(:,newVarIndex));
if rNew<=0                                                                  % check optimality
    [optimum,utility] = computeOptimum(modifiedA, b, modifiedc, base);
else
    baseSuboptimal();optimum = nan;utility = nan;
end
end

function [optimum, utility] = newCondition(modifiedA, modifiedb, c, base)
% Function that implements the routine for a new condition.
%
% Input:
%   modifiedA - nConditions+1 x nVariables
%       The matrix (normal form) with the appended new
%       condition coefficients.
%       Note, that the equations must be of "<="-type.
%   modifiedb - nConditions+1 x 1
%       The modified limits vector of the boundary conditions
%       (usually right side) with the appended new condition limit.
%   c - 1 x nVariables
%       The initial coeffient vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   base - 1 x nConditions
%       The optimal base of the initial problem.
% Output:
%   optimum - nVariables x 1
%       The optimal solution determinable by the sensitivity analysis. 
%       Note, that for invalid problems, 'nan' is returned.
%   utility - scalar
%       Value of the utility function at optimum.
base = [base size(modifiedA,2)];
if invertibleAndValid(modifiedA(:,base),modifiedb)                         % check validity
    [optimum,utility] = computeOptimum(modifiedA, modifiedb, c, base);
else
    solutionInvalid();optimum = nan;utility = nan;
end
end

function [optimum, utility] = computeOptimum(A, b, c, base)
% Function to compute the optimum w.r.t. the actual variables and utility.
c = c(1:size(A,2)-length(base));                                           % relevant part of cost vector
parameters = A(:,base)\b;                                                  % parameters optimum
optimum = zeros(length(c),1);
for i = 1:length(base)
    if base(i)<=length(optimum)
        optimum(base(i)) = parameters(i);                                  % get relevant dimensions
    end
end
utility = c*optimum;
end
    
function complement = getComplement(base, numberOfColumns)
% Compute the complement of the base
allColumns = 1:numberOfColumns;
complement = setxor(allColumns,base);
end
function base = correctBase(numberOfColumns, base)
% Correct the base indices in case a new variable has been introduced
initNumberOfVariables = numberOfColumns-length(base);
for i = 1:length(base)
    if base(i)>initNumberOfVariables
        base(i) = base(i)+1;
    end
end
end
function [standart, cost] = normal2standard(normal, cost)
% Covert normal to standard form.
numberOfConditions = size(normal,1);
standart = [normal eye(numberOfConditions)];
if nargin==2
    cost = [cost zeros(1,numberOfConditions)];
end
end

function indication = invertibleAndValid(Aj, b)
% Check if a Aj\b is computable and valid.
if det(Aj)~=0
    if min(inv(Aj)*b)>=0
        indication = 1;
        return;
    end
end
indication = 0;
solutionInvalid();
end

function baseSuboptimal()
fprintf('Error: Base suboptimal.\n');
end

function solutionInvalid()
fprintf('Error: Solution defined by base invalid.\n');
end

        
function inputError()
fprintf('Error: Please check input.\n');
end
