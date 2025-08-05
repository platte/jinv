classdef dsOperator
    % classdef dsOperator
    %
    % Authors:
    %   (c) Matthias Chung (matthias.chung@emory.edu) and Rosemary Renaut in October 2021
    %    Rosemary Renaut September 2023 
    % MATLAB Version: 9.11.0.1769968 (R2021b)
    %
    % Description:
    %   dOperator provides a class of identity, finite-difference and Laplace
    %   operators for dimension 1,2, and 3.
    %   The matrices are square 
    %
    %  D = dsOperator(type, dimension,weight)
    %
    % Properties:
    %     type       - type of matrix 'identity', 'finite difference', 'laplace', 'graph trend filter', 'diagonal', 'matrix'
    %     dimension  - provides dimension of object for matrix [m, n, o] with object dimension m x n x o
    %     weight     - allows weighting of the derivatives in each
    %     dimension
    %     sizes      - provides the dimension of the matrix D
    %     transposed - flag if operator is transposed or not
    %
    % Example:
    %     D = dsOperator('laplace',[60])
    %     t = linspace(0,2*pi,size(D,2))';
    %     x = cumsum(cumsum(sin(t)))
    %     y = D*x
    %     plot(t,x), hold on, plot(t(2:end-1),600*y)
    %
    %     D = dsOperator('finite difference',[60, 60, 60])
    %
    % References:
    %     Future versions:  use different distances in each dimension 
    %     e.g., h = [1/2 1/4 0.3452] - for the weight scaling of operator 
    %     e.g., or  mu*D (easiy to include just add at the very end a mu*(D*x) or mu*(D'*x))

    properties
        type
        dimension
        sizes
        transposed
        weight
    end

    methods

        %% initialize D operator
        function D = dsOperator(type, dimension, weight)
            if nargin < 1, version = '2025/06/09.1'; fprintf('  Current dsOperator version: %s\n', version); return, end
            D.type       = type;              % provides operator type (identity, finite difference, laplace)
            D.dimension  = dimension;         % provide dimension of operator [n1 n2 n3] refers to a 3D operator D takes inputs of length n1*n2*n3
            switch D.type                     % actual operator size of D (mxn) where n= n1*n2*n3 and m is output dimension of D, must be determined by
                case 'identity'
                    D.sizes = prod(D.dimension)*[1,1];
                case 'finite difference' %square
                    switch length(D.dimension)
                             case 1
                            D.sizes = [D.dimension, D.dimension];
                        case 2
                            D.sizes = [2*prod(D.dimension), prod(D.dimension)];
                        case 3
                            D.sizes = [3*prod(D.dimension), prod(D.dimension)];
                    %        case 1
                    %         D.sizes = [D.dimension-1, D.dimension];
                    %     case 2
                    %         D.sizes = [2*prod(D.dimension)-sum(D.dimension), prod(D.dimension)];
                    %     case 3
                    %         D.sizes = [3*prod(D.dimension)-prod(D.dimension([1,2]))-prod(D.dimension([2,3]))-prod(D.dimension([1,3])), prod(D.dimension)];
                    end
                case 'laplace' %not changed from original
                    switch length(D.dimension)
                        case 1
                            D.sizes = [D.dimension-2, D.dimension];
                        case 2
                            D.sizes = [2*prod(D.dimension)-2*sum(D.dimension), prod(D.dimension)];
                        case 3
                            D.sizes = [3*prod(D.dimension)-2*prod(D.dimension([1,2]))-2*prod(D.dimension([2,3]))-2*prod(D.dimension([1,3])), prod(D.dimension)];
                    end
                case 'diagonal'
                    if nargin < 3, error ('diagonal requires a vector ''weight''.'), end
                    D.dimension = length(weight);
                    D.weight = weight;
                    D.sizes = [D.dimension, D.dimension];
                case 'matrix'
                    if nargin < 3, error ('matrix requires a matrix ''weight''.'), end
                    D.dimension = size(weight);
                    D.weight = weight;
                    D.sizes = size(weight);
                case 'graph trend filter'
                    if nargin < 3, error ('Graph trend filter requires Adjacency Matrix ''weight''.'), end
                    D.weight = weight; % adjacency matrix
                    D.sizes = [D.dimension^2, D.dimension];
                otherwise
                    error('No valid operator type selected')
            end
            D.transposed = false;
        end

        %% transpose method
        function D = ctranspose(D)
            D.transposed = not(D.transposed);
            D.sizes = flip(D.sizes);
        end

        %% size method
        function s = size(D,dim)
            if nargin < 2
                s = D.sizes;
            else
                s = D.sizes(dim);
            end
        end

        %% mtimes method
        function x = mtimes(D, x)
            % types: identy, finite difference, and laplace
            % possible dimensions 1, 2, and 3

            if D.sizes(2)~=length(x), error('size mismatch'); end

            switch D.type
                case 'identity'
                    % nothing to do
                case 'finite difference'
                    x = applyFiniteDifferenceMatrix(D,x);
                case 'laplace'
                    x =  applyLaplaceMatrix(D,x);
                case 'diagonal'
                    x =  applyDiagonal(D,x);
                case 'matrix'
                    x =  applyMatrix(D,x);
                case 'graph trend filter'
                    x =  applyGraphTrendFilter(D,x);
            end

        end

    end
end

%%
function x = applyFiniteDifferenceMatrix(D,x)

if ~D.transposed
    switch length(D.dimension)
        case 1
            x  = diff(x);
        case 2
            x = reshape(x,D.dimension(1),D.dimension(2));
            x = [reshape(diff(x,1,1), D.dimension(2)*(D.dimension(1)-1),1); ...
                reshape(diff(x,1,2), D.dimension(1)*(D.dimension(2)-1),1)];
        case 3
            x = reshape(x,D.dimension(1),D.dimension(2),D.dimension(3));
            y1 =zeros(D.dimension(1)+1,D.dimension(2),D.dimension(3));
            y2 =zeros(D.dimension(1),D.dimension(2)+1,D.dimension(3));
            y3 =zeros(D.dimension(1),D.dimension(2),D.dimension(3)+1);
            %y =zeros(D.dimension(1)+1,D.dimension(2)+1,D.dimension(3)+1);
            y1(1:D.dimension(1), 1:D.dimension(2), 1:D.dimension(3))=x;
            y2(1:D.dimension(1), 1:D.dimension(2), 1:D.dimension(3))=x;
            y3(1:D.dimension(1), 1:D.dimension(2), 1:D.dimension(3))=x;
          %   [size(y1),size(y2),size(y3)]
            y1=diff(y1,1,1);y2=diff(y2,1,2);y3=diff(y3,1,3);
         %   [size(y1),size(y2),size(y3)]
            x=[y1(:);y2(:);y3(:)];
           % x = [reshape(diff(x,1,1), D.dimension(2)*D.dimension(3)*(D.dimension(1)-1),1); ...
            %    reshape(diff(x,1,2), D.dimension(1)*D.dimension(3)*(D.dimension(2)-1),1); ...
             %   reshape(diff(x,1,3), D.dimension(1)*D.dimension(2)*(D.dimension(3)-1),1)];
    end
else % transposed case
    switch length(D.dimension)
        case 1
            x = -diff([0; x; 0]); % augment and  difference
        case 2
            z1 = zeros(1,D.dimension(2)); % augment zeros
            z2 = zeros(D.dimension(1),1); % augment zeros
            split = D.dimension(2)*(D.dimension(1)-1);
            x1 = x(1:split);
            x2 = x(split+1:end); % split
            x1 = reshape(x1, D.dimension(1)-1, D.dimension(2) );
            x2 = reshape(x2, D.dimension(1),  D.dimension(2)-1);
            x = -reshape(diff([z1;x1;z1],1,1) + diff([z2,x2,z2],1,2),prod(D.dimension),1);
        case 3
            % z1 = zeros(1,D.dimension(2),D.dimension(3)); % augment zeros
            % z2 = zeros(D.dimension(1),1,D.dimension(3)); % augment zeros
            % z3 = zeros(D.dimension(1),D.dimension(2),1); % augment zeros
            z1 = zeros(D.dimension(1)+1,D.dimension(2),D.dimension(3)); % augment zeros
            z2 = zeros(D.dimension(1),D.dimension(2)+1,D.dimension(3)); % augment zeros
            z3 = zeros(D.dimension(1),D.dimension(2),D.dimension(3)+1); % augment zeros
            %[size(z1),size(z2),size(z3)]
            % %split1 = D.dimension(2)*D.dimension(3)*(D.dimension(1)-1);
            %split2 = D.dimension(1)*D.dimension(3)*(D.dimension(2)-1) + split1;
            split1 = D.dimension(2)*D.dimension(3)*(D.dimension(1));
            split2 = D.dimension(1)*D.dimension(3)*(D.dimension(2)) + split1;
            x1 = x(1:split1);        % split
            x2 = x(split1+1:split2); % split
            x3 = x(split2+1:end);    % split
            %[size(x1),size(x2),size(x3)]
           % x1 = reshape(x1, D.dimension(1)-1, D.dimension(2),   D.dimension(3)  );
           % x2 = reshape(x2, D.dimension(1),   D.dimension(2)-1, D.dimension(3)  );
            %x3 = reshape(x3, D.dimension(1),   D.dimension(2),   D.dimension(3)-1);
            x1 = reshape(x1, D.dimension(1), D.dimension(2),   D.dimension(3)  );
            x2 = reshape(x2, D.dimension(1),   D.dimension(2), D.dimension(3)  );
            x3 = reshape(x3, D.dimension(1),   D.dimension(2),   D.dimension(3));
           % [size(x1),size(x2),size(x3)]
            z1(1:D.dimension(1),1:D.dimension(2),1:D.dimension(3) )=x1;
            z2(1:D.dimension(1),1:D.dimension(2),1:D.dimension(3) )=x2;
            z3(1:D.dimension(1),1:D.dimension(2),1:D.dimension(3) )=x3;
           % x1 = cat(1,cat(1,z1,x1),z1); % concartinate
           % x2 = cat(2,cat(2,z2,x2),z2); % concartinate
           % x3 = cat(3,cat(3,z3,x3),z3); % concartinate
           % [size(z1),size(z2),size(z3)]
            y1=diff(z1,1,1) ;
            y2=diff(z2,1,2);
            y3=diff(z3,1,3);
            %[size(y1),size(y2),size(y3)]
            x=-reshape(y1+y2+y3,prod(D.dimension),1);
            %x = -reshape(diff(x1,1,1) ...
            %    + diff(x2,1,2) ...
            %    + diff(x3,1,3),prod(D.dimension),1);
    end
end

end

%%
function x = applyLaplaceMatrix(D,x)

if ~D.transposed
    switch length(D.dimension)
        case 1
            x  = diff(x,2);
        case 2
            x = reshape(x,D.dimension(1),D.dimension(2));
            x = [reshape(diff(x,2,1), D.dimension(2)*(D.dimension(1)-2),1); ...
                reshape(diff(x,2,2), D.dimension(1)*(D.dimension(2)-2),1)];
        case 3
            x = reshape(x,D.dimension(1),D.dimension(2),D.dimension(3));
            x = [reshape(diff(x,2,1), D.dimension(2)*D.dimension(3)*(D.dimension(1)-2),1); ...
                reshape(diff(x,2,2), D.dimension(1)*D.dimension(3)*(D.dimension(2)-2),1); ...
                reshape(diff(x,2,3), D.dimension(1)*D.dimension(2)*(D.dimension(3)-2),1)];
    end
else
    switch length(D.dimension)
        case 1
            x = -diff([0; 0; x; 0; 0],2); % augment and  difference
        case 2
            z1 = zeros(2,D.dimension(2)); % augment zeros
            z2 = zeros(D.dimension(1),2); % augment zeros
            split = D.dimension(2)*(D.dimension(1)-2);
            x1 = x(1:split);
            x2 = x(split+1:end); % split
            x1 = reshape(x1, D.dimension(1)-2, D.dimension(2) );
            x2 = reshape(x2, D.dimension(1),  D.dimension(2)-2);
            x = -reshape(diff([z1;x1;z1],2,1) + diff([z2,x2,z2],2,2),prod(D.dimension),1);
        case 3
            z1 = zeros(2,D.dimension(2),D.dimension(3)); % augment zeros
            z2 = zeros(D.dimension(1),2,D.dimension(3)); % augment zeros
            z3 = zeros(D.dimension(1),D.dimension(2),2); % augment zeros
            split1 = D.dimension(2)*D.dimension(3)*(D.dimension(1)-2);
            split2 = D.dimension(1)*D.dimension(3)*(D.dimension(2)-2) + split1;
            x1 = x(1:split1);        % split
            x2 = x(split1+1:split2); % split
            x3 = x(split2+1:end);    % split
            x1 = reshape(x1, D.dimension(1)-2, D.dimension(2),   D.dimension(3)  );
            x2 = reshape(x2, D.dimension(1),   D.dimension(2)-2, D.dimension(3)  );
            x3 = reshape(x3, D.dimension(1),   D.dimension(2),   D.dimension(3)-2);
            x1 = cat(1,cat(1,z1,x1),z1); % concartinate
            x2 = cat(2,cat(2,z2,x2),z2); % concartinate
            x3 = cat(3,cat(3,z3,x3),z3); % concartinate
            x = -reshape(diff(x1,2,1) ...
                + diff(x2,2,2) ...
                + diff(x3,2,3),prod(D.dimension),1);
    end
end
end

%%
function x = applyDiagonal(D,x)
x = D.weight.*x;
end

%%
function x = applyMatrix(D,x)
if ~D.transposed
    x = D.weight*x;
else
    x = D.weight'*x;
end
end
%%
function x =  applyGraphTrendFilter(D,x)
if ~D.transposed
    x = x.*D.weight - D.weight.*x';
    x = x(:);
else
    x = reshape(x,D.dimension*[1,1]);
    x = diag(x*D.weight' - D.weight'*x);
end
end
