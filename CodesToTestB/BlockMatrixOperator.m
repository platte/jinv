classdef BlockMatrixOperator
    % classdef Original name WGJOperator 
    %
    % Authors:
    %   (a) Rosemary Renaut in October 2022 for BTTB and BlockM
    %
    % MATLAB Version: 24.1.0.2537033 (R2024a)  
    %
    % Description:
    %  BlockMatrixOperator provides a class to implement block matrix
    %  This version concatenates matrices
    %
    %  BlockM = BlockMatrixOperator(n,blocksize,nb,blockmatrix)
    %
    % Input
    %     nb                  - number of matrices in block matrix
    %     blockmatrix    - the block matrices eg [G1;G2; ... Gnb], which may be operators
    %     scaling           - length nb for scaling the matrices by a scalar
    % Properties:
    %     sizes      - provides the dimension of the matrix G
    %     transposed - flag if operator is transposed or not
   %      nb                  - number of matrices in block matrix
    %     blockmatrix    - the block matrices eg [G1;G2; ... Gnb], which may be operators
    %     blocksizes    - the sizes of block matrices eg [G1;G2; ... Gnb], which may be operators
    %
    % Example:
    %     BlockM = BlockMatrixOperator(n,blocksize,nb,blockmatrix)
    % 
    %
    % References:
    %     Tutorial Code for BTTB

    properties
        sizes
        transposed
        nb
        blockmatrix
        blocksizes
        scaling
    end

    methods

        %% initialize BlockM operator
        function BlockM = BlockMatrixOperator(nb,blockmatrix,scaling)
            if nargin < 1, version = '2025/03/08.1';
                fprintf('  Current BlockMOperator version: %s\n', version);
                return,
            end
            for j=1:nb
                blocksizes(j,:)=blockmatrix{j}.sizes;
            end
            BlockM.sizes                = sum(blocksizes,1);
            BlockM.transposed      = false;
            BlockM.nb                    = nb;
            BlockM.blocksizes       = blocksizes;
            BlockM.blockmatrix      = blockmatrix;
            BlockM.scaling             = scaling;
        end

        %% transpose method
        function BlockM = ctranspose(BlockM)
            BlockM.transposed = not(BlockM.transposed);
            BlockM.sizes = flip(BlockM.sizes);
        end

        %% size method
        function s = size(BlockM,dim)
            if nargin < 2
                s = BlockM.sizes;
            else
                s = BlockM.sizes(dim);
            end
        end

        %% mtimes method
        function y = mtimes(BlockM, x)
            if BlockM.sizes(2)~=length(x(:)),
                error('size mismatch');
            end
            j=1;j1=1;
            aa=BlockM.blockmatrix{j};
            if ~BlockM.transposed 
                j2=BlockM.blocksizes(j,2);
                y=(aa*x(j1:j2,1))*BlockM.scaling(j);
                for j=2:BlockM.nb
                    aa=BlockM.blockmatrix{j};
                    j1=j2+1;j2=BlockM.blocksizes(j,2)+j2;
                    y=[y;(aa*x(j1:j2,1))*BlockM.scaling(j)];
                end
            else
                j2=BlockM.blocksizes(j,1);
                y=(aa'*x(j1:j2,1))*BlockM.scaling(j);
                for j=2:BlockM.nb
                    aa=BlockM.blockmatrix{j};
                    j1=j2+1;j2=BlockM.blocksizes(j,1)+j2;
                    y=[y;(aa'*x(j1:j2,1))*BlockM.scaling(j);];
                end
            end
        end
    end
end
%%