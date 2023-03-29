classdef Data < handle

    properties
        SubjectID
        Name
        Gender
        Weight

        %Landmark
        REFM
        REFL
        RMM
        RLM

        %Joint center
        RHipCenter
        RKneeCenter
        RAnkleCenter

        %HipFitting
        hipFitErrorMean
        hipFitErrorSD




        %外部尺寸
        Height%身高
        RShankL%小腿外侧长
        RShankM%小腿内侧长
        RKneeW%膝盖宽
        RAnkleW%膝盖宽

        %body frame
        REFMInit
        REFLInit
        RMMInit
        RLMInit
        RKneeCenterInit
        RAnkleCenterInit



        %JCS坐标系
        RFemoralRotm
        RKneeRotm%大腿坐标系下
        RTibiaRotm
        RFemoralEul
        RTibiaEul
        RKneeEul

        RHipQuat
        RKneeQuat
        RFemoralQuat
        RTibiaQuat


        RFemoralClusterEul
        RTibiaClusterEul
        RFemoralClusterRotm
        RTibiaClusterRotm
        RFemoralClusterPos
        RTibiaClusterPos

        %膝关节
        RKneeAoR_F
        RKneeAoR_T
        RKneeAoR_G

        RKneeAoR_F_Selected
        RKneeAoR_T_Selected

        MeanRKneeAoR_F
        MeanRKneeAoR_T

        SDRKneeAoR_F
        SDRKneeAoR_T

        AngleAoR_FYZ
        AngleAoR_FXZ
        AngleAoR_TYZ
        AngleAoR_TXZ


        frameForUnify
        ReferenceFrameInit
        RotmF
        RotmT
        RotmT2

        ErrorFemoralRotm
        ErrorTibiaRotm

        ROMF
        ROMT

        temp



    end

    methods

        function obj = Data()
        end

        function clusterEul2Rotm(obj)
            RFemoralClusterRotmTemp =  eul2rotm(obj.RFemoralClusterEul);
            RTibiaClusterRotmTemp =  eul2rotm(obj.RTibiaClusterEul);

            obj.RFemoralClusterRotm =  cell([size(RFemoralClusterRotmTemp,3),1]);
            obj.RTibiaClusterRotm =  cell([size(RFemoralClusterRotmTemp,3),1]);

            for i = 1:size(RFemoralClusterRotmTemp,3)
                obj.RFemoralClusterRotm{i,1} = RFemoralClusterRotmTemp(:,:,i);
                obj.RTibiaClusterRotm{i,1} = RTibiaClusterRotmTemp(:,:,i);
            end
        end

        function fillMissingPoint(obj)
            obj.REFM = fillmissing(obj.REFM,"spline");
            obj.REFL = fillmissing(obj.REFL,"spline");
            obj.RMM = fillmissing(obj.RMM,"spline");
            obj.RLM = fillmissing(obj.RLM,"spline");
            obj.RFemoralClusterEul= fillmissing(obj.RFemoralClusterEul,"spline");
            obj.RTibiaClusterEul = fillmissing(obj.RTibiaClusterEul,"spline");

        end

        function [fitErrorMean,fitErrorSD] = setHipAsOrigin(obj)
            data=[obj.REFM;obj.REFL];
            f=@(p,data)(data(:,1)-p(1)).^2+(data(:,2)-p(2)).^2+(data(:,3)-p(3)).^2-p(4)^2;
            [p,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(data,zeros(size(data,1),1),f,[0 0 0 1]');%拟合的参数
            hipPosition = p(1:3);
            obj.REFM = obj.REFM-hipPosition';
            obj.REFL = obj.REFL-hipPosition';
            obj.RMM = obj.RMM-hipPosition';
            obj.RLM = obj.RLM-hipPosition';
            obj.RKneeCenter = (obj.REFM+obj.REFL)/2;
            obj.RAnkleCenter = (obj.RMM+obj.RLM)/2;
            obj.RHipCenter = zeros(size(obj.REFM,1),3);

            if size(obj.RFemoralClusterPos) ~= 0

                obj.RFemoralClusterPos = obj.RFemoralClusterPos-hipPosition';
                obj.RTibiaClusterPos = obj.RTibiaClusterPos-hipPosition';
            end


            r = abs(p(4));
            %             [X,Y,Z] = sphere;
            %             X2 = X * r;
            %             Y2 = Y * r;
            %             Z2 = Z * r;
            %             surf(X2,Y2,Z2)

            data=[obj.REFM;obj.REFL];
            error = sqrt(sum(data.^2')')-r;
            obj.hipFitErrorMean = mean(error);
            obj.hipFitErrorSD = std(error);
            fitErrorMean = obj.hipFitErrorMean;
            fitErrorSD = obj.hipFitErrorSD;

        end
        function setHipAsOriginXsens(obj)
            obj.REFM = obj.REFM-obj.RHipCenter;
            obj.REFL = obj.REFL-obj.RHipCenter;
            obj.RMM = obj.RMM-obj.RHipCenter;
            obj.RLM = obj.RLM-obj.RHipCenter;
            obj.RKneeCenter = (obj.REFM+obj.REFL)/2;
            obj.RAnkleCenter = (obj.RMM+obj.RLM)/2;
            obj.RHipCenter = zeros(size(obj.REFM,1),3);
        end

        function rotateReferenceFrame(obj,Rotm)

            if size(obj.REFM,1) ~= 0
                obj.REFM = (Rotm*obj.REFM')';
                obj.REFL = (Rotm*obj.REFL')';
                obj.RMM = (Rotm*obj.RMM')';
                obj.RLM = (Rotm*obj.RLM')';
                obj.RKneeCenter = (Rotm*obj.RKneeCenter')';
                obj.RAnkleCenter = (Rotm*obj.RAnkleCenter')';
            end

            if size(obj.RFemoralClusterPos) ~= 0
                obj.RFemoralClusterPos = (Rotm*obj.RFemoralClusterPos')';
                obj.RTibiaClusterPos = (Rotm*obj.RTibiaClusterPos')';
                for i = 1:size(obj.REFM,1)
                    obj.RFemoralClusterRotm{i} = Rotm*obj.RFemoralClusterRotm{i};
                    obj.RFemoralClusterRotm{i} = Rotm*obj.RFemoralClusterRotm{i};
                end
            end

            if size(obj.RFemoralRotm,1) ~= 0
                for i = 1:size(obj.RFemoralRotm,1)
                    obj.RFemoralRotm{i} = Rotm*obj.RFemoralRotm{i};
                end
            end

            if size(obj.RTibiaRotm,1)~= 0
                for i = 1:size(obj.RTibiaRotm,1)
                    obj.RTibiaRotm{i} = Rotm*obj.RTibiaRotm{i};
                end
            end
        end

        function getFemoralOrientation(obj)
            for i = 1:size(obj.RKneeCenter,1)
                %x: The line perpendicular to plane defined bythe origin and the two FEs
                % both y- and z-axis, pointing anteriorly
                xAxis = cross(obj.REFL(i,:)-obj.RHipCenter(i,:),obj.REFM(i,:)-obj.RHipCenter(i,:));
                xAxis = xAxis/norm(xAxis);
                %y: The line joining the midpoint between the medial and lateral FEs and the origin, and pointing cranially.
                yAxis = -obj.RKneeCenter(i,:)+obj.RHipCenter(i,:);
                yAxis = yAxis/norm(yAxis);
                %z: The line perpendicular to both y- and x-axis, pointing to the right
                zAxis = cross(xAxis,yAxis);
                zAxis = zAxis/norm(zAxis);

                obj.RFemoralRotm{i,1} = [xAxis',yAxis',zAxis'];
                obj.RFemoralEul(i,:) = rotm2eul(obj.RFemoralRotm{i,1});
            end
        end

        function getTibiaOrientationXsens(obj)
            %Since the knee rotation axis could not be calculated through dynamic data, we adopted the knee definition of Xsens MVN.
            %Y: Ankle center to knee center
            %Z: EFM to EFL, pointing right
            for i = 1:size(obj.RKneeCenter,1)
                yAxis = obj.RKneeCenter(i,:)-obj.RAnkleCenter(i,:);
                yAxis = yAxis/norm(yAxis);
                zAxis = (obj.REFL-obj.REFM);
                zAxis = zAxis/norm(zAxis);
                xAxis = cross(yAxis,zAxis);
                xAxis = xAxis/norm(xAxis);
                obj.RTibiaRotm{i,1} = [xAxis',yAxis',zAxis'];
                obj.RTibiaEul(i,:) = rotm2eul(obj.RTibiaRotm{i,1});
            end
        end

        function getKneeAoR_F_old(obj,startFrame,frameStep,endFrame,clusterFlag)
            threshold = 0.5;
            j = 1;
            for i = startFrame:frameStep:endFrame
                k = i;
                l = i+frameStep;
                if clusterFlag
                    RKneeRotm_G = obj.RTibiaClusterRotm{l}/obj.RTibiaClusterRotm{k};
                else
                    RKneeRotm_G = obj.RTibiaRotm{l}/obj.RTibiaRotm{k};
                end

                obj.RKneeAoR_G(j,1) = i;
                obj.RKneeAoR_G(j,2:5) = rotm2axang(RKneeRotm_G);

                obj.RKneeAoR_F(j,1) = i;
                obj.RKneeAoR_F(j,2:4) = (obj.RFemoralRotm{k}'*obj.RKneeAoR_G(j,2:4)')';

                if obj.RKneeAoR_F(j,4)<0
                    obj.RKneeAoR_F(j,2:4) = -obj.RKneeAoR_F(j,2:4);
                end
                obj.RKneeAoR_F(j,5) = obj.RKneeAoR_G(j,5);
                j = j+1;
            end
            obj.RKneeAoR_G(any(isnan(obj.RKneeAoR_G), 2),:) = [];
            obj.RKneeAoR_F(any(isnan(obj.RKneeAoR_F), 2),:) = [];

            maxAngle = max(obj.RKneeAoR_G(:,5));
            k = 1;
            for i= 1:size(obj.RKneeAoR_F,1)
                if obj.RKneeAoR_F(i,5)>threshold*maxAngle
                    obj.RKneeAoR_F_Selected(k,:) = obj.RKneeAoR_F(i,:);
                    k = k+1;
                end
            end
            OA = mean(obj.RKneeAoR_F_Selected(:,2:4).*obj.RKneeAoR_F_Selected(:,5));
            obj.MeanRKneeAoR_F = OA/norm(OA);
            re = obj.RKneeAoR_F_Selected(:,2:4)-obj.MeanRKneeAoR_F;
            for i = 1:size(re)
                a(i) = asin(re(i)/2)*2 ./pi*180;
            end
            obj.SDRKneeAoR_F = std(a);
        end

        function [calculationFlag,window] = getAoRCalculationWindow(obj)
            maxPos = max(obj.REFM(:,3));
            minPos = min(obj.REFM(:,3));
            trajectoryRange = maxPos-minPos;
            minThreshold = (minPos+0.15*trajectoryRange).*[1,1];
            maxThreshold = (maxPos-0.15*trajectoryRange).*[1,1];

            %             t = [1,3000];
            %             figure
            %             hold on
            %             plot(t,minThreshold)
            %             plot(t,maxThreshold)
            %             plot(obj.REFM(:,3))


            calculationFlag = zeros(size(obj.REFM,1),1);
            windowFlag = 0;
            window = [];
            windowCount = 1;
            for j = 1:size(obj.REFM,1)
                if obj.REFM(j,3)>minThreshold(1)&&obj.REFM(j,3)<maxThreshold(1)
                    calculationFlag(j) = 1;
                end
            end

            for j = 1:size(obj.REFM,1)-1
                if calculationFlag(j) == 0 && calculationFlag(j+1) == 1
                    window(1,windowCount) = j;
                    windowCount = windowCount+1;
                end
            end

            windowCount = 1;
            for j = 1:size(obj.REFM,1)-1
                if calculationFlag(j) == 1 && calculationFlag(j+1) == 0
                    window(2,windowCount) = j;
                    windowCount = windowCount+1;
                end
            end
            window(3,:) = window(2,:)-window(1,:);
        end

        function getKneeAoR_F(obj,window,stepCount, clusterFlag)
            j = 1;
            for i = 1:size(window,2)
                start = window(1,i);
                step = round(window(3,i)/stepCount);
                if step<10
                    continue
                end

                for k = start:step:start+step*stepCount
                    if k+step>size(obj.REFM,1)
                        continue
                    end

                    if clusterFlag
                        RKneeRotm_F1 = obj.RFemoralRotm{k}'*obj.RTibiaClusterRotm{k};
                        RKneeRotm_F2 = obj.RFemoralRotm{k+step}'*obj.RTibiaClusterRotm{k+step};
                        RKneeRotm_F = RKneeRotm_F2/RKneeRotm_F1;
                    else
                        %                         RKneeRotm_F = obj.RTibiaRotm{k+step}/obj.RFemoralRotm{k};
                        RKneeRotm_F1 = obj.RFemoralRotm{k}'*obj.RTibiaRotm{k};
                        RKneeRotm_F2 = obj.RFemoralRotm{k+step}'*obj.RTibiaRotm{k+step};
                        RKneeRotm_F = RKneeRotm_F2/RKneeRotm_F1;
                    end

                    obj.RKneeAoR_F(j,1) = k;
                    obj.RKneeAoR_F(j,2:5) = rotm2axang(RKneeRotm_F);
                    obj.RKneeAoR_G(j,1) = k;
                    obj.RKneeAoR_G(j,2:4) = obj.RFemoralRotm{k}*obj.RKneeAoR_F(j,2:4)';

                    if obj.RKneeAoR_F(j,4)<0
                        obj.RKneeAoR_F(j,2:4) = -obj.RKneeAoR_F(j,2:4);
                    end

                    j = j+1;
                end
            end

            angleThreshold = 0.5;
            zThreshold = 0.6;

            maxAngle = max(obj.RKneeAoR_F(:,5));
            k = 1;
            for i= 1:size(obj.RKneeAoR_F,1)
                if obj.RKneeAoR_F(i,5)>angleThreshold*maxAngle&&obj.RKneeAoR_F(i,4)>zThreshold
                    obj.RKneeAoR_F_Selected(k,:) = obj.RKneeAoR_F(i,:);
                    k = k+1;
                end
            end

            OA = mean(obj.RKneeAoR_F_Selected(:,2:4).*obj.RKneeAoR_F_Selected(:,5));
            obj.MeanRKneeAoR_F = OA/norm(OA);
            re = obj.RKneeAoR_F_Selected(:,2:4)-obj.MeanRKneeAoR_F;
            for i = 1:size(re)
                a(i) = asin(re(i)/2)*2 ./pi*180;
            end
            for i = 1:size(obj.RKneeAoR_F_Selected,1)
                a(i) = acos(dot(obj.RKneeAoR_F_Selected(i,2:4),obj.MeanRKneeAoR_F)/ ...
                    (norm(obj.RKneeAoR_F_Selected(i,2:4))*norm(obj.MeanRKneeAoR_F)))./pi*180;
            end
            obj.SDRKneeAoR_F = std(a);
            obj.AngleAoR_FYZ = atan(obj.MeanRKneeAoR_F(:,2)./obj.MeanRKneeAoR_F(:,3))./pi*180;
            obj.AngleAoR_FXZ = atan(obj.MeanRKneeAoR_F(:,1)./obj.MeanRKneeAoR_F(:,3))./pi*180;

        end

        function getTibiaCluster(obj)
            for i = 1:size(obj.RKneeCenter,1)
                xAxis = cross(obj.RLM(i,:)-obj.RKneeCenter(i,:),obj.RMM(i,:)-obj.RKneeCenter(i,:));
                xAxis = xAxis/norm(xAxis);
                %y: The line joining the midpoint between the medial and lateral FEs and the origin, and pointing cranially.
                yAxis = -obj.RAnkleCenter(i,:)+obj.RKneeCenter(i,:);
                yAxis = yAxis/norm(yAxis);
                %z: The line perpendicular to both y- and x-axis, pointing to the right
                zAxis = cross(xAxis,yAxis);
                zAxis = zAxis/norm(zAxis);
                obj.RTibiaClusterRotm{i,1} = [xAxis',yAxis',zAxis'];
            end
        end


        function getTibiaOrientation(obj)
            %y:line between knee center and the ankle joint centre
            YAxis = obj.RKneeCenter-obj.RAnkleCenter;
            for i = 1:size(obj.RMM,1)
                yAxis = YAxis(i,:)/norm(YAxis(i,:));
                Axis_G = obj.RFemoralRotm{i}*obj.MeanRKneeAoR_F';  %knee AoR in reference frame
                %x:line perpendicular to both yAxis and Axis_G
                xAxis = cross(yAxis,Axis_G);
                xAxis = xAxis/norm(xAxis);

                zAxis = cross(xAxis,yAxis);
                zAxis = zAxis/norm(zAxis);

                obj.RTibiaRotm{i,1} = [xAxis',yAxis',zAxis'];
                obj.RTibiaEul(i,:) = rotm2eul(obj.RTibiaRotm{i,1});
                obj.RKneeAoR_T(i,:) = obj.RTibiaRotm{i,1}'*Axis_G;
            end
            obj.RKneeAoR_T(any(isnan(obj.RKneeAoR_T), 2),:) = [];
            obj.MeanRKneeAoR_T = mean(obj.RKneeAoR_T);
        end

        function getKneeAoR_T(obj)
            for i = 1:size(obj.RFemoralRotm,1)
                Axis_G = obj.RFemoralRotm{i}*obj.MeanRKneeAoR_F';  %全局膝关节轴
                obj.RKneeAoR_T(i,:) = obj.RTibiaRotm{i,1}'*Axis_G;
            end
            obj.RKneeAoR_T(any(isnan(obj.RKneeAoR_T), 2),:) = [];
            obj.MeanRKneeAoR_T = mean(obj.RKneeAoR_T);

            for i = 1:size(obj.RKneeAoR_T,1)
                a(i) = acos(dot(obj.RKneeAoR_T(i,1:3),obj.MeanRKneeAoR_T)/ ...
                    (norm(obj.RKneeAoR_T(i,1:3))*norm(obj.MeanRKneeAoR_T)))./pi*180;
            end
            obj.SDRKneeAoR_T = std(a);


            obj.AngleAoR_TYZ = atan(obj.MeanRKneeAoR_T(:,2)./obj.MeanRKneeAoR_T(:,3))./pi*180;
            obj.AngleAoR_TXZ = atan(obj.MeanRKneeAoR_T(:,1)./obj.MeanRKneeAoR_T(:,3))./pi*180;


        end

        function getMarkPositionInSegment(obj)
            obj.REFMInit = (obj.RFemoralRotm{1}'*(obj.REFM(1,:)-obj.RHipCenter(1,:))')';
            obj.REFLInit = (obj.RFemoralRotm{1}'*(obj.REFL(1,:)-obj.RHipCenter(1,:))')';
            obj.RMMInit = (obj.RTibiaRotm{1}'*(obj.RMM(1,:)-obj.RKneeCenter(1,:))')';
            obj.RLMInit = (obj.RTibiaRotm{1}'*(obj.RLM(1,:)-obj.RKneeCenter(1,:))')';
            obj.RKneeCenterInit = (obj.REFMInit + obj.REFLInit)/2;
        end

        function getBodyDimension(obj)
            obj.RShankM = norm(obj.REFMInit-obj.RMMInit-obj.RKneeCenterInit);
            obj.RShankL = norm(obj.REFLInit-obj.RLMInit-obj.RKneeCenterInit);
            obj.RKneeW = norm(obj.REFMInit-obj.REFLInit);
            obj.RAnkleW = norm(obj.RMMInit-obj.RLMInit);
            obj.RKneeCenterInit = norm(obj.REFMInit+obj.REFLInit)/2;
            obj.RAnkleCenterInit = norm(obj.RMMInit+obj.RLMInit)/2;
        end


        function rotateFemoralFrame(obj,Rotm)
            for i = 1:size(obj.RFemoralRotm,1)
                obj.RFemoralRotm{i} = obj.RFemoralRotm{i}*Rotm;
            end
        end

        function rotateTibiaFrameInFemoralFrame(obj,Rotm)
            for i = 1:size(obj.RFemoralRotm,1)
                obj.RTibiaRotm{i} = Rotm*obj.RTibiaRotm{i};
            end
        end

        function rotateTibiaFrame(obj,Rotm)
            for i = 1:size(obj.RFemoralRotm,1)
                obj.RTibiaRotm{i} = obj.RTibiaRotm{i}*Rotm;
            end
        end


        function getRKneeRotm(obj)
            for i = 1:size(obj.RMM,1)
            end
        end





        function getMarkPosition(obj)
            if obj.Gender
                load("Class\Data\MarkTableM.mat")
                for i = 1:20
                    if obj.SubjectID == MarkTableM.SubjectID(i)

                        obj.REFMInit = [MarkTableM.REFMx(i),MarkTableM.REFMy(i),MarkTableM.REFMz(i)]';
                        obj.REFLInit = [MarkTableM.REFLx(i),MarkTableM.REFLy(i),MarkTableM.REFLz(i)]';
                        obj.RMMInit = [MarkTableM.RMMx(i),MarkTableM.RMMy(i),MarkTableM.RMMz(i)]';
                        obj.RLMInit = [MarkTableM.RLMx(i),MarkTableM.RLMy(i),MarkTableM.RLMz(i)]';
                        return
                    end
                end
            else
                load("Class\Data\MarkTableF.mat")
                for i = 1:18
                    if obj.SubjectID == MarkTableF.SubjectID(i)
                        obj.REFMInit = [MarkTableF.REFMx(i),MarkTableF.REFMy(i),MarkTableF.REFMz(i)]';
                        obj.REFLInit = [MarkTableF.REFLx(i),MarkTableF.REFLy(i),MarkTableF.REFLz(i)]';
                        obj.RMMInit = [MarkTableF.RMMx(i),MarkTableF.RMMy(i),MarkTableF.RMMz(i)]';
                        obj.RLMInit = [MarkTableF.RLMx(i),MarkTableF.RLMy(i),MarkTableF.RLMz(i)]';
                        return
                    end
                end
            end
        end


        function useOMCMarkPosition(obj,omcData)
            obj.REFMInit = omcData.REFMInit;
            obj.REFLInit = omcData.REFLInit;
            obj.RMMInit = omcData.RMMInit;
            obj.RLMInit = omcData.RLMInit;
        end


        function getMarkerPositionOMC(obj)
            if obj.Gender
                load("Class\Data\MarkTableM.mat")
                for i = 1:20
                    if obj.SubjectID == MarkTableM.SubjectID(i)
                        obj.REFMInit = [MarkTableM.REFMX(i),MarkTableM.REFMY(i),MarkTableM.REFMZ(i)]';
                        obj.REFLInit = [MarkTableM.REFLX(i),MarkTableM.REFLY(i),MarkTableM.REFLZ(i)]';
                        obj.RMMInit = [MarkTableM.RMMX(i),MarkTableM.RMMY(i),MarkTableM.RMMZ(i)]';
                        obj.RLMInit = [MarkTableM.RLMX(i),MarkTableM.RLMY(i),MarkTableM.RLMZ(i)]';
                        return
                    end
                end
            else
                load("Class\Data\MarkTableF.mat")
                for i = 1:18
                    if obj.SubjectID == MarkTableF.SubjectID(i)
                        obj.REFMInit = [MarkTableF.REFMX(i),MarkTableF.REFMY(i),MarkTableF.REFMZ(i)]';
                        obj.REFLInit = [MarkTableF.REFLX(i),MarkTableF.REFLY(i),MarkTableF.REFLZ(i)]';
                        obj.RMMInit = [MarkTableF.RMMX(i),MarkTableF.RMMY(i),MarkTableF.RMMZ(i)]';
                        obj.RLMInit = [MarkTableF.RLMX(i),MarkTableF.RLMY(i),MarkTableF.RLMZ(i)]';
                        return
                    end
                end
            end



        end

        function getMarkTrajectories(obj)
            for i = 1:size(obj.RFemoralRotm,1)
                obj.REFL(i,:) = (obj.RFemoralRotm{i,1}*obj.REFLInit)';
                obj.REFM(i,:) = (obj.RFemoralRotm{i,1}*obj.REFMInit)';
                obj.RMM(i,:) = (obj.RTibiaRotm{i,1}*obj.RMMInit)';
                obj.RLM(i,:) = (obj.RTibiaRotm{i,1}*obj.RLMInit)';
            end
            obj.RKneeCenter = (obj.REFL+obj.REFM)/2;
            obj.RMM =  obj.RMM+obj.RKneeCenter;
            obj.RLM = obj.RLM+obj.RKneeCenter;
            obj.RAnkleCenter = (obj.RMM+obj.RLM)/2;
        end

        function [RotmF,RotmT] = getCorrectionRotm(obj,MeanRKneeAoR_F,MeanRKneeAoR_T)
            angle = acos( ...
                dot(obj.MeanRKneeAoR_F,MeanRKneeAoR_F)/ ...
                (norm(obj.MeanRKneeAoR_F)*norm(MeanRKneeAoR_F)));
            axis = cross(obj.MeanRKneeAoR_F,MeanRKneeAoR_F);
            axis = axis/norm(axis);
            RotmF = axang2rotm([axis,angle]);

            angle = acos( ...
                dot(obj.MeanRKneeAoR_T,MeanRKneeAoR_T)/ ...
                (norm(obj.MeanRKneeAoR_T)*norm(MeanRKneeAoR_T)));
            axis = cross(obj.MeanRKneeAoR_T,MeanRKneeAoR_T);
            axis = axis/norm(axis);
            RotmT = axang2rotm([axis,angle]);
        end

        function getInitFrame(obj)
            %             data = [obj.RKneeCenter;obj.RAnkleCenter];
            obj.RKneeCenter = (obj.REFL+obj.REFM)/2;
            data = [obj.RKneeCenter];
            f=@(p,data)(data(:,1)*p(1))+(data(:,2)*p(2))+(data(:,3)*p(3));
            [p,~,~,~,~,~] = nlinfit(data,zeros(size(data,1),1),f,[1 0 0]');%拟合的参数

            z = obj.RFemoralRotm{1}(:,3);
            if dot(z,p)<0
                ZAxis = -p/norm(p);
            else
                ZAxis = p/norm(p);
            end

            YAxis = -obj.RKneeCenter(1,:);
            YAxis = YAxis/norm(YAxis);

            XAxis = cross(YAxis,ZAxis);
            XAxis = XAxis/norm(XAxis);

            YAxis = cross(ZAxis,XAxis);
            YAxis = YAxis/norm(YAxis);
            obj.ReferenceFrameInit = [XAxis',YAxis',ZAxis];

        end

        function getInitFrameByAxis(obj)
            j = 1;
            for i = 1:100:1000
                k = i;
                l = i+100;
                RHipRotm_G = obj.RFemoralRotm{l}/obj.RFemoralRotm{k};
                RHipAoR_G(j,1) = i;
                RHipAoR_G(j,2:5) = rotm2axang(RHipRotm_G);
                j = j+1;
            end

            for i = 1:size(RHipAoR_G,1)
                frame = RHipAoR_G(i,1);
                flag = dot(RHipAoR_G(i,2:4),obj.RFemoralRotm{frame}(:,3));
                if flag<0
                    RHipAoR_G(i,2:4) = -RHipAoR_G(i,2:4);
                end
            end
            obj.ZAxis = mean(RHipAoR_G(:,2:4));
        end

        function getFrameForUnifyByPeakPoint(obj,bottom,top)
            ZAxis = cross(obj.RKneeCenter(bottom,:),obj.RKneeCenter(top,:));
            YAxis = obj.RKneeCenter(bottom,:);
            ZAxis = ZAxis/norm(ZAxis);
            YAxis = YAxis/norm(YAxis);
            XAxis = cross(YAxis,ZAxis);
            XAxis = XAxis/norm(XAxis);
            obj.frameForUnify = [XAxis',YAxis',ZAxis'];
        end

        function  getFemoralFrameCorrRotm(obj,MeanAoRF)
            angle = acos( ...
                dot(MeanAoRF,obj.MeanRKneeAoR_F)/ ...
                (norm(MeanAoRF)*norm(obj.MeanRKneeAoR_F)));
            Axis = cross(MeanAoRF,obj.MeanRKneeAoR_F);
            Axis = Axis/norm(Axis);
            a = angle./pi*180;
            axang = [Axis,angle];
            obj.RotmF = axang2rotm(axang);
            %             if abs(angle)<30/180*pi;
            %                 obj.RotmF = eye(3);
            %             end


        end

        function  getTibiaFrameCorrRotm(obj,MeanAoRT)
            angle = acos( ...
                dot(MeanAoRT,obj.MeanRKneeAoR_T)/ ...
                (norm(MeanAoRT)*norm(obj.MeanRKneeAoR_T)));
            Axis = cross(MeanAoRT,obj.MeanRKneeAoR_T);
            Axis = Axis/norm(Axis);
            axang = [Axis,angle];
            a = angle/pi*180;
            obj.RotmT = axang2rotm(axang);
            %             if abs(angle)<30/180*pi;
            %                 obj.RotmT = eye(3);
            %             end
        end

        function getTibiaFrameCorrRotm2(obj,omcInitRotm)

            obj.RotmT2 = obj.RTibiaRotm{1}'*omcInitRotm;
            %             obj.RotmT2 = obj.RTibiaRotm{1}^-1*omcInitRotm;
            %             obj.RotmT2 = omcInitRotm/obj.RTibiaRotm{1};
            %             obj.RotmT2 = obj.RTibiaRotm{1}/omcInitRotm;
        end

        function getEul(obj)
            %             disp("debug")
            %             RotInit = [-1 0 0;0 0 1;0 1 0];
            %             RotInit = [0 0 -1;-1 0 0;0 1 0];

            for i = 1:size(obj.RFemoralRotm,1)
                obj.RFemoralEul(i,:) = rotm2eul(obj.RFemoralRotm{i,1});
                obj.RTibiaEul(i,:) = rotm2eul(obj.RTibiaRotm{i,1});
                %         RKneeEul
%                 RomF = RotInit'*obj.RFemoralRotm{i,1};
%                 RFemoralEul_init(i,:) = rotm2eul(RomF);
%                 RomT = RotInit'*obj.RTibiaRotm{i,1};
%                 RTibiaEul_init(i,:) = rotm2eul(RomT);
            end

            %             obj.ROMF(1) = max(RFemoralEul_init(:,1))-min(RFemoralEul_init(:,1));
            %             obj.ROMF(2) = max(RFemoralEul_init(:,2))-min(RFemoralEul_init(:,2));
            %             obj.ROMF(3) = max(RFemoralEul_init(:,3))-min(RFemoralEul_init(:,3));
            %             obj.ROMT(1) = max(RTibiaEul_init(:,1))-min(RTibiaEul_init(:,1));
            %             obj.ROMT(2) = max(RTibiaEul_init(:,2))-min(RTibiaEul_init(:,2));
            %             obj.ROMT(3) = max(RTibiaEul_init(:,3))-min(RTibiaEul_init(:,3));


            %             obj.ROMF(1) = max(abs(obj.RFemoralEul(:,1)));
            %             obj.ROMF(2) = max(abs(obj.RFemoralEul(:,2)));
            %             obj.ROMF(3) = max(abs(obj.RFemoralEul(:,3)));
            %             obj.ROMT(1) = max(abs(obj.RTibiaEul(:,1)));
            %             obj.ROMT(2) = max(abs(obj.RTibiaEul(:,2)));
            %             obj.ROMT(3) = max(abs(obj.RTibiaEul(:,3)));

            max(obj.REFM(3));
            topFrame = find(obj.REFM(:,3)==max(obj.REFM(:,3)));
            a = [0,0,-1];
            b = obj.RFemoralRotm{topFrame(1),1}(:,2);
            obj.ROMF = 180-acos(dot(a,b)/(norm(a)*norm(b)))./pi*180;
%             a = [0,0,-1];
%             b = imcData.RTibiaRotm{topFrame(1),1}(:,2);
%             obj.ROMT = acos(dot(a,b)/(norm(a)*norm(b)))./pi*180;
            obj.ROMT = 90;




        end

        function xsensResample(obj,originSampleRate)
            obj.REFM = resample(obj.REFM,100,originSampleRate);
            obj.REFL = resample(obj.REFL,100,originSampleRate);
            obj.RMM = resample(obj.RMM,100,originSampleRate);
            obj.RLM = resample(obj.RLM,100,originSampleRate);
            obj.RHipCenter = resample(obj.RHipCenter,100,originSampleRate);
            %             obj.RFemoralQuat = resample(obj.RFemoralQuat,100,originSampleRate);
            %             obj.RTibiaQuat = resample(obj.RTibiaQuat,100,originSampleRate);

        end

        function getXsensLandmarks(obj,xsensMarkers)
            obj.REFM = xsensMarkers.pRightKneeMedEpicondyle;
            obj.REFL = xsensMarkers.pRightKneeLatEpicondyle;
            obj.RMM = xsensMarkers.pRightMedMalleolus;
            obj.RLM = xsensMarkers.pRightLatMalleolus;
        end

        function getXsensQuat(obj,xsensQuat)
            obj.RFemoralQuat = xsensQuat{:,62:65};
            obj.RTibiaQuat = xsensQuat{:,66:69};
        end

        function getXsensRotm(obj)

            RFemoralRotmTemp = quat2rotm(obj.RFemoralQuat);
            RTibiaRotmTemp = quat2rotm(obj.RTibiaQuat);
            obj.RFemoralRotm =  cell([size(RFemoralRotmTemp,3),1]);
            obj.RTibiaRotm =  cell([size(RFemoralRotmTemp,3),1]);

            for i = 1:size(RFemoralRotmTemp,3)
                obj.RFemoralRotm{i,1} = RFemoralRotmTemp(:,:,i);
                obj.RTibiaRotm{i,1} = RTibiaRotmTemp(:,:,i);
            end

        end





        function getOMCLandmarks(obj,omcMarkers)
            if size(fieldnames(omcMarkers),1)==4
                obj.REFM = omcMarkers.REFM;
                obj.REFL = omcMarkers.REFL;
                obj.RMM = omcMarkers.RMM;
                obj.RLM = omcMarkers.RLM;
            else
                obj.REFM = omcMarkers.iW3_Digitized1;
                obj.REFL = omcMarkers.iW3_Digitized2;
                obj.RMM = omcMarkers.iW2_Digitized1;
                obj.RLM = omcMarkers.iW2_Digitized2;
            end
            obj.REFM(obj.REFM==0) = NaN;
            obj.REFL(obj.REFL==0) = NaN;
            obj.RMM(obj.RMM==0) = NaN;
            obj.RLM(obj.RLM==0) = NaN;
        end


        function delectC3DMissingPoints(obj)
            obj.REFM(obj.REFM==0) = NaN;
            obj.REFL(obj.REFL==0) = NaN;
            obj.RMM(obj.RMM==0) = NaN;
            obj.RLM(obj.RLM==0) = NaN;
            %             obj.REFM = fillmissing(obj.REFM,"spline");
            %             obj.REFL = fillmissing(obj.REFL,"spline");
            %             obj.RMM = fillmissing(obj.RMM,"spline");
            %             obj.RLM = fillmissing(obj.RLM,"spline");
        end

        function selectRange(obj,startFrame,endFrame)
            obj.REFM = obj.REFM(startFrame:endFrame,:);
            obj.REFL = obj.REFL(startFrame:endFrame,:);
            obj.RMM = obj.RMM(startFrame:endFrame,:);
            obj.RLM = obj.RLM(startFrame:endFrame,:);
            if size(obj.RHipCenter) ~= 0
                obj.RHipCenter = obj.RHipCenter(startFrame:endFrame,:);
                obj.RKneeCenter = obj.RKneeCenter(startFrame:endFrame,:);
                obj.RAnkleCenter = obj.RAnkleCenter(startFrame:endFrame,:);
            end

            if size(obj.RFemoralRotm,1)~=0
                obj.RFemoralRotm = obj.RFemoralRotm(startFrame:endFrame,:);
                obj.RTibiaRotm = obj.RTibiaRotm(startFrame:endFrame,:);
            end
        end
    end
end