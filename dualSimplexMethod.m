% Last modified: 28.06.2017


function [optimum, utility, numberOfIterations] = dualSimplexMethod(conditionCoefficients, limits, utilityCoefficients, info, startingBase)
% This is a function that implements the  dual simplex method.
% Note that any solutions where xi < 0 are disregarded.

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
%       The optimal solution determinable by the simplex method. Note, that
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
if (nargin==5)                                                             % if given, utilize starting base
    tableau = changeCompleteBase(tableau, startingBase);
    base = startingBase;
else
    base = defaultBase;
end
    
    
numberOfIterations = 0;
while true
    numberOfIterations = numberOfIterations+1;
    pivot = getNewPivot(tableau,base,info);                                % get new pivot element
    [optimum,utility] = getParameterValues(tableau,base);                  % compute current parameter values and utility
    if info,optionalPrint(tableau,utility,pivot);end                       % optional print
    if isnan(pivot),return;end                                             % if no valid pivot was found terminate
    base(pivot(1))= pivot(2);                                              % update base
    tableau = changeBase(tableau, pivot);                                  % update tableau by doing gauss elimination
end
end

function tableau = changeCompleteBase(tableau, base)
% Change all elements of the base at once.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%   base - 1 x nConditions
%       The current base vector, containing the valid columns in the order
%       they make up the identity matrix.
%
% Output:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The tableau corresponding to the new base.
for i = 1:length(base)
    tableau = changeBase(tableau, [i base(i)]);                      % iterate over all elements of the base
end
end

function optionalPrint(tableau,utility,pivot)
% Print the optional iteration output.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%   utility - scalar
%       Value of the utility function at current solution.
%   pivot - 1 x 2
%       The pivot element specified by [rowIndex columnIndex]
%       in the tableau.
tableau
if ~isnan(pivot)
    fprintf('Utility: %f \n',utility);
end
end

function [parameterValues, utility] = getParameterValues(tableau,base)
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
for column = 1:numberOfParameters                                          % match the parameters to base
    row = find(base==column);
    if isempty(row)
        parameterValues(column) = 0;
    else
        parameterValues(column) = tableau(row,end);
    end
end
utility = -tableau(end,end);
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

function pivot = getNewPivot(tableau, base, info)
% A new pivot element is determined or the optimization is terminated.
%
% Input:
%   tableau - nConditions+1 x nParameters+nConditions+2
%       The current simplex tableau.
%   base - 1 x nConditions
%       The current base vector, containing the valid columns in the order
%       they make up the identity matrix.
%   info - scalar (bool)
%       Enables/Disables informative output.
%
% Ouput:
%   pivot - 1 x 2
%       The pivot element specified by [rowIndex columnIndex]
%       in the tableau.

pivot = 0;
while pivot==0                                                             % while there are possible pivot elements
    [candidateRowValue,candidateRowIndex] = min(tableau(1:end-1,end));
    if (candidateRowValue<0)                                               % check whether valid
        rowIndex = candidateRowIndex;
        candidateRow = tableau(rowIndex,1:end-2);                          % set candidate row
        candidateColumnIndex = minOfTheNegatives(candidateRow,tableau(end,1:end-2),base);  % set candidate column index
        if isnan(candidateColumnIndex)                                     % if not set to dummy value (note, that this is a internal copy of the tableau)
            tableau(rowIndex,end)=Inf;
        else
            columnIndex = candidateColumnIndex;
            pivot = [rowIndex columnIndex];                                % else set
        end
    else
        if info,disp('No further improvements possible:');end              % if no improvement is possible, terminate (further error detection is disregarded for readability)
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

function index = minOfTheNegatives(row,costRow,base)
% Searches for columnindex of the next pivot element.
% Input:
%   row - 1 x nParameters+nConditions
%       The candidate row.
%   costRow - 1 x nParameters+nConditions
%       The cost vector of the standard form.
%
%   base - 1 x nConditions
%       The current base vector, containing the valid columns in the order
%       they make up the identity matrix.
%
% Output:
%   index - scalar
%       The index corresponding to the selected pivot's column.
%       Defaults to "nan" if only invalid elements are contained.
costRow(base) = nan;                                                       
indexVector = 1:length(row);
indexVector = indexVector(row<0);                                          % get indices of negative elements
if isempty(indexVector),index = nan;return;end                             % if only invalid elements..
[~,I] = min(abs(costRow(indexVector)./row(indexVector)));                  % get relative index of selected element
index = indexVector(I);                                                    % absolute index
end
        


        

        
