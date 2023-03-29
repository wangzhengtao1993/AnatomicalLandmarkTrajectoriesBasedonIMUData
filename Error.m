classdef Error < handle
    %UNTITLED2 此处提供此类的摘要
    %   此处提供详细说明

    properties
        iREFM
        iREFL
        iRMM
        iRLM

        oREFM
        oREFL
        oRMM
        oRLM

        REFM
        REFL
        RMM
        RLM

        rmsREFM
        rmsREFL
        rmsRMM
        rmsRLM

        rREFM
        rREFL
        rRMM
        rRLM

        ICCREFM
        ICCREFL
        ICCRMM
        ICCRLM
        icc

        pearson
        RMS

        rmsFemoralAngle
        rmsTibiaAngle

        rmsFemoralYAxis
        rmsTibiaYAxis



    end

    methods
        function obj = Error(imcData,omcData)
            obj.iREFM = imcData.REFM;
            obj.iREFL = imcData.REFL;
            obj.iRMM = imcData.RMM;
            obj.iRLM = imcData.RLM;

            obj.oREFM = omcData.REFM;
            obj.oREFL = omcData.REFL;
            obj.oRMM = omcData.RMM;
            obj.oRLM = omcData.RLM;
            obj.getError;
            obj.getRMS;
            obj.getCorr;
            obj.getICC;

        end

        function getError(obj)
            obj.REFM = obj.iREFM-obj.oREFM;
            obj.REFL = obj.iREFL-obj.oREFL;
            obj.RMM = obj.iRMM-obj.oRMM;
            obj.RLM = obj.iRLM-obj.oRLM;
            for i = 1:size(obj.REFM,1)
                obj.REFM(i,4) = norm(obj.REFM(i,1:3));
                obj.REFL(i,4) = norm(obj.REFL(i,1:3));
                obj.RMM(i,4) = norm(obj.RMM(i,1:3));
                obj.RLM(i,4) = norm(obj.RLM(i,1:3));
            end
        end

        function getCorr(obj)
            obj.iREFM(any(isnan(obj.oREFM), 2),:) = [];
            obj.iREFL(any(isnan(obj.oREFL), 2),:) = [];
            obj.iRMM(any(isnan(obj.oRMM), 2),:) = [];
            obj.iRLM(any(isnan(obj.oRLM), 2),:) = [];

            obj.oREFM(any(isnan(obj.oREFM), 2),:) = [];
            obj.oREFL(any(isnan(obj.oREFL), 2),:) = [];
            obj.oRMM(any(isnan(obj.oRMM), 2),:) = [];
            obj.oRLM(any(isnan(obj.oRLM), 2),:) = [];


            for i = 1:3
                obj.rREFM(i) = corr(obj.iREFM(:,i),obj.oREFM(:,i));
                obj.rREFL(i) = corr(obj.iREFL(:,i),obj.oREFL(:,i));
                obj.rRMM(i) = corr(obj.iRMM(:,i),obj.oRMM(:,i));
                obj.rRLM(i) = corr(obj.iRLM(:,i),obj.oRLM(:,i));
            end

            obj.pearson = [obj.rREFM,obj.rREFL,obj.rRMM,obj.rRLM];
            %             obj.pearson = obj.rREFM(i)
        end

        function getICC(obj)
            obj.iREFM(any(isnan(obj.oREFM), 2),:) = [];
            obj.iREFL(any(isnan(obj.oREFL), 2),:) = [];
            obj.iRMM(any(isnan(obj.oRMM), 2),:) = [];
            obj.iRLM(any(isnan(obj.oRLM), 2),:) = [];

            obj.oREFM(any(isnan(obj.oREFM), 2),:) = [];
            obj.oREFL(any(isnan(obj.oREFL), 2),:) = [];
            obj.oRMM(any(isnan(obj.oRMM), 2),:) = [];
            obj.oRLM(any(isnan(obj.oRLM), 2),:) = [];
            type = 'C-1';


            for i = 1:3
                [r, LB, UB, F, df1, df2, p] = ICC([obj.iREFM(:,i),obj.oREFM(:,i)], type);
                obj.ICCREFM(i) = r;
                [r, LB, UB, F, df1, df2, p] = ICC([obj.iREFL(:,i),obj.oREFL(:,i)], type);
                obj.ICCREFL(i) = r;
                [r, LB, UB, F, df1, df2, p] = ICC([obj.iRMM(:,i),obj.oRMM(:,i)], type);
                obj.ICCRMM(i) = r;
                [r, LB, UB, F, df1, df2, p] = ICC([obj.iRLM(:,i),obj.oRLM(:,i)], type);
                obj.ICCRLM(i) = r;
            end

            obj.icc = [obj.rREFM,obj.rREFL,obj.rRMM,obj.rRLM];
            %             obj.pearson = obj.rREFM(i)
        end


        function getRMS(obj)
            obj.REFM(any(isnan(obj.REFM), 2),:) = [];
            obj.REFL(any(isnan(obj.REFL), 2),:) = [];
            obj.RMM(any(isnan(obj.RMM), 2),:) = [];
            obj.RLM(any(isnan(obj.RLM), 2),:) = [];
            
            obj.rmsREFM = rms(obj.REFM);
            obj.rmsREFL = rms(obj.REFL);
            obj.rmsRMM = rms(obj.RMM);
            obj.rmsRLM = rms(obj.RLM);
            obj.RMS = [obj.rmsREFM(4),obj.rmsREFL(4),obj.rmsRMM(4),obj.rmsRLM(4)];
        end

        function rmsTable = getRMSTable(obj,RMSArray)
            MeanRMS = mean(RMSArray);
            SDRMS = std(RMSArray);

            for i = 2:5
                RMSTable{i-1} = [num2str(MeanRMS(i),'%.1f'),char(177),num2str(SDRMS(i),'%.1f')];
            end
            rmsTable = cell2table(RMSTable);
            rmsTable.Properties.VariableNames = {'REFM','REFL','RMM','RLM'};
        end

        function iccTable = getICCTable(obj,ICCArray,ICCArrayM,ICCArrayF)
            ICCArray(:,1) = [];
            ICCArrayM(:,1) = [];
            ICCArrayF(:,1) = [];

            MAX = max(ICCArrayM,[],'all');
            MIN = min(ICCArrayM,[],'all');
            M = [num2str(MAX,'%.3f'),'~',num2str(MIN,'%.3f')];
            MAX = max(ICCArrayF,[],'all');
            MIN = min(ICCArrayF,[],'all');
            F = [num2str(MAX,'%.3f'),'~',num2str(MIN,'%.3f')];
            MAX = max(ICCArray,[],'all');
            MIN = min(ICCArray,[],'all');
            T = [num2str(MAX,'%.3f'),'~',num2str(MIN,'%.3f')];
            iccArray = [M;F;T]                     

            iccTable = cell2table(iccArray);
%             rmsTable.Properties.VariableNames = {'REFM','REFL','RMM','RLM'};
        end




        function rTable = getRTable(obj,CorrArray)
            MeanRMS = mean(CorrArray);
            SDRMS = std(CorrArray);
            for i = 2:13
                rTable{i-1} = [num2str(MeanRMS(i),'%.3f'),char(177),num2str(SDRMS(i),'%.3f')];
            end
        end

        function getAngleRMS(obj,omcData,imcData)
            k = 1;
            for j = 1:size(imcData.RFemoralRotm,1)
                if isnan(omcData.RFemoralRotm{j,1}(1,1))
                    continue
                end
                errorFemoralRotm{k,1} = imcData.RFemoralRotm{j,1}/omcData.RFemoralRotm{j,1};
                if isnan(omcData.RTibiaRotm{j,1}(1,1))
                    continue
                end
                errorTibiaRotm{k,1} = imcData.RTibiaRotm{j,1}/omcData.RTibiaRotm{j,1};
                k = k+1;
            end

            for j = 1:size(errorFemoralRotm,1)
                errorFemoralEul(j,:) = rotm2eul(errorFemoralRotm{j})/pi*180;
            end

            for j = 1:size(errorTibiaRotm,1)

                errorTibiaEul(j,:) = rotm2eul(errorTibiaRotm{j})/pi*180;
            end

            obj.rmsFemoralAngle = rms(errorFemoralEul);
            obj.rmsTibiaAngle = rms(errorTibiaEul);

        end
        
        function getYAxisAngleRMS(obj,omcData,imcData)
            k = 1;
            for j = 1:size(imcData.RFemoralRotm,1)
                if isnan(omcData.RFemoralRotm{j,1}(1,1))
                    continue
                end
                imcFemoralYAxis = imcData.RFemoralRotm{j,1}(:,2);
                omcFemoralYAxis = omcData.RFemoralRotm{j,1}(:,2);
                errorFemoralYAxis(k) = acos(dot(imcFemoralYAxis,omcFemoralYAxis))./pi*180;

                if isnan(omcData.RTibiaRotm{j,1}(1,1))
                    continue
                end
                imcTibiaYAxis = imcData.RTibiaRotm{j,1}(:,2);
                omcTibiaYAxis = omcData.RTibiaRotm{j,1}(:,2);
                errorTibiaYAxis(k) = acos(dot(imcTibiaYAxis,omcTibiaYAxis))./pi*180;
                k = k+1;
            end

            obj.rmsFemoralYAxis = rms(errorFemoralYAxis);
            obj.rmsTibiaYAxis = rms(errorTibiaYAxis);
           
        end



        


    end
end