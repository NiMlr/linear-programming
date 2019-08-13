% Last modified: 14.06.2017


function [optimum, utility, numberOfIterations] = simplexMethod(conditionCoefficients, limits, utilityCoefficients, info, startingBase)
% This is a function that implements the simplex method.
% Note, that the origin must be a valid solution or a valid starting base
% needs to be given. Note that any solutions where xi < 0 are disregarded.

% Input:
%   conditionCoefficients - nConditions x nParameters
%       A matrix with the coeffitients (usually left side) of the boundary
%       condition-equations (along the lines). Note, that the equations
%       must be of "<="-type.
%   limits - nConditions x 1
%       The limits vector of the boundary conditions (usually right side).
%   utilityCoefficients - 1 x nParameters
%       This is the coeffient Vector of the linear utility-function.
%       Note, that this is a maximising algorithm.
%   info - scalar (bool)
%       Enables/Disables informative output.
%   startingBase (optional) - 1 x numberOfConditions
%       Provide a starting base to start optimization process from.
%       If this input is omitted, the origin will be used instead.
%
% Output:
%   optimum - nParameters x 1
%       The best solution determinable by the simplex method. Note, that
%       for invalid problems, the best solution found is returned.
%   utility - scalar
%       Value of the utility function at optimum.
%   numberOfIterations - scalar
%       Number of iterations till expected termination.
% ------------------------------------------------------------------------
% Example input:
%   Conditions: 6 x1 + 15x2 <= 4500                            [ 6 15;
%               4 x1 + 5 x2 <= 2000 -> conditionCoeffitients =   4  5;
%               20x1 + 10x2 <= 8000                             20 10 ]
%                               |             [ 4500;
%                                --> limits =   2000;
%                                               8000 ];
%   Utility-function: max 16x1 + 32x2  ->  utilityCoeffitients = [ 16 32 ];
%   You might wish to take a look at whats happening: info = true;
%   You could provide the starting base: startingBase = [ 3 4 5 ];
% ------------------------------------------------------------------------

[tableau, defaultBase] = setUpTableau(conditionCoefficients, limits, utilityCoefficients); % set up the tableau
if (nargin==5)                                                             % check whether starting base is given and it is valid
    tableau = checkValidity(tableau,startingBase);                         % change tableau accordingly
    if isnan(tableau),return;end
    base = startingBase;
else,base = defaultBase;                                                   % otherwise use default base
end
numberOfIterations = 0;
while true
    numberOfIterations = numberOfIterations+1;
    pivot = getNewPivot(tableau);                                          % get new pivot element
    [optimum,utility] = getParameterValues(tableau,base);                  % compute current parameter values and utility
    if info,fprintf('Current solution: %s \n Utility: %f \n',mat2str(optimum),utility);end   % optional print
    if isnan(pivot),return;end                                             % if no valid pivot was found terminate
    base(pivot(1))= pivot(2);                                              % update base
    tableau = changeBase(tableau, pivot);                                  % update tableau by doing gauss elimination
end
end


function [parameterValues, utility] = getParameterValues(tableau, base)
% Function that computes current parameter values by stepping through the
% base vector and checking the according row of the last column of the
% tableau.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%   base - 1 x nConditions
%       The current base vector, containing the valid columns in the order
%       they make up the identity matrix.
%
% Output:
%   parameterValues - nParameters x 1
%       The current best parameter values.
%   utility - scalar
%       Utility of current best parameter values.
numberOfParameters = size(tableau,2)-2-length(base);
parameterValues = zeros(numberOfParameters,1);
for column = 1:numberOfParameters
    row = find(base==column);
    if isempty(row)
        parameterValues(column) = 0;
    else
        parameterValues(column) = tableau(row,end);
    end
end
utility = -tableau(end,end);
end
    

function tableau = checkValidity(tableau, startingBase)
% Checks validity of the tableau initialized with an optional given base by
% checking for invalid values after performing gauss elimination.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The initial simplex tableau.
%   startingBase - The optional non-identity (slack variables) base.
%
% Output:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       Will be returned after performing gauss elimination and checking
%       it. If it turns out to be invalid, an error is generated and the
%       tableau is set to nan.

if length(unique(startingBase))==length(startingBase)
    for i = 1:length(startingBase)
        tableau = changeBase(tableau, [i startingBase(i)]);                % pivot over all base elements (rows specified by order)
    end
    for i = 1:size(tableau,1)-1
        if (isnan(tableau(i,end))||tableau(i,end)<0)                       % check values
            disp(tableau)
            disp('Invalid base')
            tableau = nan;return;
        end
    end
end
end

function tableau = changeBase(tableau, pivot)
% Change the base by adding a new column and omitting another.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%   pivot - 1 x 2
%       The pivot element specified by [rowIndex columnIndex]
%       in the tableau.
%
% Output:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The tableau corresponding to the new base.

pivotRow = pivot(1);
pivotColumn = pivot(2);

tableau(pivotRow,:) = tableau(pivotRow,:)/tableau(pivotRow,pivotColumn);   % gauss elimination
for i = 1:size(tableau,1)
   if i ~= pivotRow 
       tableau(i,:) = tableau(i,:)- tableau(i,pivotColumn)*tableau(pivotRow,:);
   end
end
end

function pivot = getNewPivot(tableau)
% A new pivot element is determined or the optimization is terminated.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%
% Ouput:
%   pivot - 1 x 2
%       The pivot element specified by [rowIndex columnIndex]
%       in the tableau.

pivot = 0;
while pivot==0                                                             % while there are possible pivot elements
    [candidateColumnShadowCost,candidateColumnIndex] = max(tableau(end,1:end-2));
    if (candidateColumnShadowCost>0)                                       % check whether valid
        columnIndex = candidateColumnIndex;
        column = tableau(1:end-1,columnIndex);
        limits = tableau(1:end-1,end);
        rowIndex = minPos(limits./column);
        isError = checkPosZero(column, tableau(1:end-1,end));
        if isnan(rowIndex)||isError                                        % if not set to dummy value (note, that this is a internal copy of the tableau)
            tableau(end,columnIndex)=-Inf;
        else
            pivot = [rowIndex columnIndex];                                % else set
        end
    else
        disp('No further improvements possible:')                          % if no improvement is possible, terminate (further error detection is disregarded for readability)
        pivot = nan;
    end
end

end

function [tableau, defaultBase] = setUpTableau(conditionCoefficients, limits, utilityCoefficients)
% Function that sets up the internal tableau representation.
%
% Inputs:
%   conditionCoefficients - nConditions x nParameters
%       A matrix with the coeffitients (usually left side) of the boundary
%       condition-equations (along the lines).
%   limits - nConditions x 1
%       The limits vector of the boundary conditions (usually right side).
%   utilityCoefficients - 1 x nParameters
%       This is the coeffient Vector of the linear utility-function.
%
% Outputs:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The initial simplex tableau.
%   defaultBase - 1 x nConditions
%       Vector that contains column indices of current base in order.

numberOfConditions = size(conditionCoefficients,1);
numberOfParameters = size(conditionCoefficients,2);

slackAndUtility = eye(numberOfConditions+1);                               % set slack part of tableau
tableau = [[[conditionCoefficients;utilityCoefficients] slackAndUtility] [limits;0]];       % put together tableau

defaultBase = numberOfParameters+1:numberOfParameters+numberOfConditions;  % generate default base
end

function I = minPos(vector)
% Function that determines the index of the minimum postive value in
% "vector". (Can be improved to be O(n) instead of probably ~O(nlogn)
% currently)

[vector,oriIndi] = sort(vector);
for i = 1:length(vector)
    if vector(i)>0
        I=oriIndi(i);return;
    end
end
I = nan;
end
        
function isError = checkPosZero(column, b)
% Check for "positive Zero" Error. Which would lead to negative Values.
for i = 1:length(column)
    if (column(i)>0)&&(b(i)==0)
        isError = 1;
        return;
    end
    isError = 0;
end

end
        
