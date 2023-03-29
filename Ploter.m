classdef Ploter < handle
    %UNTITLED7 此处提供此类的摘要
    %   此处提供详细说明

    properties
        SubjectID
        angle_YZ_sd
        angle_XZ_sd
        angle_YZ
        angle_XZ
        meanAxis
    end

    methods
        function obj = Ploter(SubjectID)
            obj.SubjectID = SubjectID;
        end

        function plotInit(obj)
            figure
            title(obj.SubjectID)
            set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.5);
            set(gcf,'unit','centimeters','position',[1,1,8,8]);set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.5);
            hold on;axis equal;box on;grid on;view([0 90])
            xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
            xlim([-800 200])
            ylim([-647 400])
            zlim([-468 600])
            view([-153 31])
            view([-130 31])
        end


        function plotInitTPose(obj)
            figure
            title(obj.SubjectID)
            set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.5);
            set(gcf,'unit','centimeters','position',[1,1,8,8]);set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.5);
            hold on;axis equal;box on;grid on;view([0 90])
            xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)');
            xlim([-327 268])
            ylim([-286 184])
            view([48 23])

        end

        function plotFemoralFrameInit(obj)
            quiver3(0,0,0,200,0,0,'g',LineWidth=1)
            quiver3(0,0,0,0,200,0,'b',LineWidth=1)
            quiver3(0,0,0,0,0,200,'r',LineWidth=1)
        end

        function plotMarkerTrajectories(obj,Data,lineSpec,startFrame,frameStep,endFrame)
            plot3(Data.REFM(startFrame:frameStep:endFrame,1),...
                Data.REFM(startFrame:frameStep:endFrame,2),...
                Data.REFM(startFrame:frameStep:endFrame,3),lineSpec)

            plot3(Data.REFL(startFrame:frameStep:endFrame,1),...
                Data.REFL(startFrame:frameStep:endFrame,2),...
                Data.REFL(startFrame:frameStep:endFrame,3),lineSpec)

            plot3(Data.RMM(startFrame:frameStep:endFrame,1),...
                Data.RMM(startFrame:frameStep:endFrame,2),...
                Data.RMM(startFrame:frameStep:endFrame,3),lineSpec)

            plot3(Data.RLM(startFrame:frameStep:endFrame,1),...
                Data.RLM(startFrame:frameStep:endFrame,2),...
                Data.RLM(startFrame:frameStep:endFrame,3),lineSpec)
        end

        function plotFemoralFrame(obj,Data,startFrame,frameStep,endFrame)
            for i = startFrame:frameStep:endFrame
                obj.plotFrame([0,0,0],Data.RFemoralRotm{i},100)
                quiver3(0,0,0,-Data.RFemoralRotm{i}(1,2),...
                    -Data.RFemoralRotm{i}(2,2),...
                    -Data.RFemoralRotm{i}(3,2),400,'g')
            end
        end

        function plotFemoralFrameClusterFrame(obj,Data,startFrame,frameStep,endFrame)
            for i = startFrame:frameStep:endFrame
                obj.plotFrame(Data.RFemoralClusterPos(i,:),Data.RFemoralClusterRotm{i},100)
                obj.plotFrame([-600,300,0],Data.RFemoralClusterRotm{i},100)

            end
        end

        function plotTibiaFrameClusterFrame(obj,Data,startFrame,frameStep,endFrame)
            for i = startFrame:frameStep:endFrame
                obj.plotFrame(Data.RTibiaClusterPos(i,:),Data.RTibiaClusterRotm{i},100)
%                 obj.plotFrame([-600,0,0],Data.RTibiaClusterRotm{i},100)
            end
        end



        function plotFemoralFrameOMC(obj,Data,startFrame,frameStep,endFrame)
            for i = startFrame:frameStep:endFrame

                obj.plotFrame([0,0,0],Data.RFemoralRotm{i},100)

                quiver3(0,0,0,-Data.RFemoralRotm{i}(1,2),...
                    -Data.RFemoralRotm{i}(2,2),...
                    -Data.RFemoralRotm{i}(3,2),400,'g')
            end
        end

        function plotTibiaFrame(obj,Data,startFrame,frameStep,endFrame)
            for i = startFrame:frameStep:endFrame
                obj.plotFrame(Data.RKneeCenter(i,:),Data.RTibiaRotm{i},100)

                quiver3(Data.RKneeCenter(i,1),Data.RKneeCenter(i,2),Data.RKneeCenter(i,3), ...
                    -Data.RTibiaRotm{i}(1,2),...
                    -Data.RTibiaRotm{i}(2,2),...
                    -Data.RTibiaRotm{i}(3,2),400,'g')

%                 obj.plotFrame([0,300,0],Data.RTibiaRotm{i},100)
            end
        end

        function plotFrame(obj,P,R,scale)
            quiver3(P(1),P(2),P(3),R(1,1),R(2,1),R(3,1),scale,'r')
            quiver3(P(1),P(2),P(3),R(1,2),R(2,2),R(3,2),scale,'g')
            quiver3(P(1),P(2),P(3),R(1,3),R(2,3),R(3,3),scale,'b')
        end

        function plotKneeAoR_T(obj,Data,scale)
            for i = 1:size(Data.RKneeAoR_T_Selected,1)

                quiver3(0,0,0,Data.RKneeAoR_T_Selected(i,2), ...
                    Data.RKneeAoR_T_Selected(i,3), ...
                    Data.RKneeAoR_T_Selected(i,4),scale,'b')

            end

        end

        function plotKneeAoR(obj,Axes)
            figure
            hold on;axis equal;box on;grid on;
            set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.5);
            set(gcf,'unit','centimeters','position',[1,1,16,12]);
            xlim([-200,200]);ylim([-200,200]);zlim([-60,200])
            xlabel('X (mm)');ylabel('Y (mm)');zlabel('Z (mm)')
            view([-35 31])
            MeanAxis = mean(Axes);
            for i = 1:size(Axes,1)
                quiver3(0,0,0,Axes(i,1),Axes(i,2),Axes(i,3),200,'b')
                text(200*Axes(i,1),200*Axes(i,2),200*Axes(i,3),num2str(Axes(i,4)))
            end
            x = MeanAxis(1);
            y = MeanAxis(2);
            z = MeanAxis(3);
            obj.meanAxis = MeanAxis;



            angle_YZ = atan(Axes(:,2)./Axes(:,3))./pi*180;
            angle_XZ = atan(Axes(:,1)./Axes(:,3))./pi*180;

            obj.angle_YZ_sd = std(angle_YZ);
            obj.angle_XZ_sd = std(angle_XZ);
            obj.angle_YZ = mean(angle_YZ);
            obj.angle_XZ = mean(angle_XZ);
        end

        function Error = plotErrorVector(obj,imcData,omcData)

            Error.REFM = imcData.REFM-omcData.REFM;
            Error.REFL = imcData.REFL-omcData.REFL;
            Error.RMM = imcData.RMM-omcData.RMM;
            Error.RLM = imcData.RLM-omcData.RLM;

            for i = 1:100:1000
                quiver3(imcData.REFM(i,1),imcData.REFM(i,2),imcData.REFM(i,3),...
                    Error.REFM(i,1),Error.REFM(i,2),Error.REFM(i,3),'off','m')
                quiver3(imcData.REFL(i,1),imcData.REFL(i,2),imcData.REFL(i,3),...
                    Error.REFL(i,1),Error.REFL(i,2),Error.REFL(i,3),'off','m')
                quiver3(imcData.RMM(i,1),imcData.RMM(i,2),imcData.RMM(i,3),...
                    Error.RMM(i,1),Error.RMM(i,2),Error.RMM(i,3),'off','m')
                quiver3(imcData.RLM(i,1),imcData.RLM(i,2),imcData.RLM(i,3),...
                    Error.RLM(i,1),Error.RLM(i,2),Error.RLM(i,3),'off','m')
            end

            
        end

        function plotValue(obj,imcData,omcData)
            figure
            subplot(2,2,1)
            hold on
            plot(imcData.REFM(:,1),'r')
            plot(omcData.REFM(:,1),'r--')

        end


    end
end