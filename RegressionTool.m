classdef RegressionTool < handle

    properties
        Property1
    end

    methods
        function obj = RegressionTool()

        end

        function t_RegressEqu = getRegressEquOld(obj,x,y,rowName)
            RegressEqu = zeros(1,5);
            mdl = fitlm(x,y);
            e = mdl.Residuals{:,1};
            if size(x,2) ==1
                RegressEqu(1,1:2) = mdl.Coefficients{:,1}';
            else
                RegressEqu(1,1:3) = mdl.Coefficients{:,1}';
            end

            RegressEqu(1,4) = mdl.Rsquared.Ordinary;
            RegressEqu(1,5) = mdl.ModelFitVsNullModel.Pvalue;
            RegressEqu(1,6) = std(e);

            t_RegressEqu = array2table(RegressEqu);
            t_RegressEqu.Properties.VariableNames = {'b0','b1','b2','RSqure','p-Value','SD'};
            t_RegressEqu.Properties.RowNames = {rowName};
        end

        function RegressEqu = getRegressEqu(obj,x,y)
            RegressEqu = zeros(1,5);
            mdl = fitlm(x,y);
            e = mdl.Residuals{:,1};
            if size(x,2) ==1
                RegressEqu(1,1:2) = mdl.Coefficients{:,1}';
            else
                RegressEqu(1,1:3) = mdl.Coefficients{:,1}';
            end

            RegressEqu(1,4) = mdl.Rsquared.Ordinary;
            RegressEqu(1,5) = mdl.ModelFitVsNullModel.Pvalue;
            RegressEqu(1,6) = std(e);
        end




        function error = getTestError_b0b1(obj,x,y,RegressEqu)
            error = y-(RegressEqu(1)+x*RegressEqu(2)); 
        end


        function error = getTestError_b0b1b2(obj,x,y,RegressEqu)
            error = y-(RegressEqu(1)+x(:,1)*RegressEqu(2)+x(:,2)*RegressEqu(3));
        end

        function errorPer = getTestErrorPer_b0b1(obj,x,y,RegressEqu)
            error = y-(RegressEqu(1)+x*RegressEqu(2));
            errorPer = abs(error./y)*100;
        end


        function errorPer = getTestErrorPer_b0b1b2(obj,x,y,RegressEqu)
            error = y-(RegressEqu(1)+x(:,1)*RegressEqu(2)+x(:,2)*RegressEqu(3));
            errorPer = abs(error./y)*100;

        end

        function [meanError,SDError] = getTestErrorMeanSD(obj,testError)       
%             if Gender == 'm'
%                 k = 20;
%             else
%                 k =18;
%             end
            k = size(testError,1)*size(testError,2);
            meanError = mean(reshape(testError,[k,size(testError,3)]));
            SDError = std(reshape(testError,[k,size(testError,3)]));
        end

        function regressTable = getRegressTable(obj,RegressEqu,RegressEqu_RowNames)
            regressTable = array2table(RegressEqu);
            regressTable.Properties.VariableNames = {'b0','b1','b2','RSqure','p-Value','SD'};
            regressTable.Properties.RowNames = RegressEqu_RowNames;

        end

        function errorTable = errorArray2Table(obj,meanError,SDError)
        
        end


    end
end