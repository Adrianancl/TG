classdef TemperatureAnalyze
%TEMPERATUREANALYZE Build an object with general informations about one
%previously determined matrix from an entire struct of these matrices.
%OBJ = TEMPERATUREANALYZE(FRAME,INDEX) returns an object with all the
%properties of this class.
    
    properties(GetAccess = 'public', SetAccess = 'private')
        CoordinateToolTip;
        TemperatureToolTip;
        RakeAngle;%Rake face slope
        ClearanceAngle;%Clearance face slope
        ShearAngle;
        FrictionAngle;
        MeanTemperatureTool;
        MaximumTemperatureTool;
        MaximumTemperatureChip;
        MaximumTemperatureCuttingZone;
        HeatCarriedAwayByChip;
        HeatFluxAwayFromToolTip;
        HeatFluxThroughWorkpiece;
        TotalPowerBalance;%Total energy produced during the process
        InternalEnergyTool;%Energy only for valid pixels
        CuttingForcePowerDirection;
        CuttingForceUncutChipThicknessDirection;        
        CuttingForceParallelToolFace;
        CuttingForceParallelShearPlane;
        CuttingForcePerpendicularShearPlane;
        CuttingForcePerpendicularToolFace;
        CoefficientFriction;
        ShearStress;
        NormalStress;
        PecletNumber;
        RatioR;
        ShearEnergyVolume;
        FrictionEnergyVolume;
        CuttingVelocity;
        UnCutChipThickness;
        ContactLength;
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        coordRF;
        coordCF;
        BW;
        lines;
        frame;%Variable of the single frame that will be analyzed
        pointCF;%auxiliar to plot the cutting edge
        pointRF;
        pointM;%Point on the bissectrix of RF and CF
        Tx;%auxiliar to plot the gradients of the frame x axis
        Ty;%auxiliar to plot the gradients of the frame y axis
        biImageTool;%Binary image of the tool shape
        biImageChip;
        biShearLine;
        xyMaxTemp;%coordinates of the point inside the chip with maximum Temperature
        lineChip;
        lineTool;%
        validTemperature;%Minimum temperature inside the valid range
        heatCapacityTool;%Heat capacity for the tool(regression of points)
        heatCapacityWP;%Heat capacity for the tool(regression of points)
        heatConductivity;%Heat conductivity for the tool
        nExcPoints;
        heatAccumulatedPerLine;
        ptosLines;
        extPtosLineChip;
        line200;%auxiliar variable to get the last valid isotherm
    end
    
    methods
    %Here are listed the functions that are used to build the object    
        function obj = TemperatureAnalyze(Frame,index)%constructor
            %Inputs -----------------------------------------
            Fp =  Frame.Fp;%Cutting force in the power direction (Newtons)
            Fq = Frame.Fq;%Passive force (Newtons)
            widthTool = Frame.w;%meters
            ap = Frame.ap;%meters
            Vp = Frame.Vp/60;%meters/second
            tuc = Frame.tuc;%meters
            clength = Frame.clength;%contact length tool - chip           
            obj.validTemperature = 200;% For any experiment
            A = 0.1;%percentage of the deformation energy that is converted in heat
            obj.heatConductivity = 75.4;%Tool
            obj.heatCapacityWP = [-4.39956806034758e-07 0.000707314520321484 ...
            -0.0488770693887544 481.214007868631]; %AISI 1045
            obj.heatCapacityTool = (10^6)*[2.50542895559373e-10 ...
            -1.99579761670655e-06 0.00274369536032376 3.09265830398264];
            %-------------------------------------------------
            obj.CuttingVelocity = Vp*60;%m/minute
            obj.UnCutChipThickness = tuc;
            obj.frame = Frame.f(index).f;
            obj.ClearanceAngle = Frame.alpha;
            obj.RakeAngle = Frame.gamma;
            if isempty(clength)
                clength = obj.contactLength();
            end
            obj.ContactLength = clength;
            if isempty(Frame.ToolTip)
                obj = obj.calculateCoordinates();
                if ~isempty(obj.coordRF) && ~isempty(obj.coordCF)
                    obj = obj.coordinateToolTip();
                else %Default conditions
                    figure
                    imagesc(Frame.f(index).f)
                    Frame.ToolTip = input('The code was not able to find the tool tip position. Provide valid coordinates [x y]: ');
                    close
                    obj.ClearanceAngle = Frames.alpha;
                    obj.RakeAngle = Frames.gamma;
                end
            end
            obj.CoordinateToolTip = Frame.ToolTip;            
            obj = obj.toolContour();
            obj = obj.findLineTool();
            obj = obj.chipContour();
            obj = obj.findLineChip();
%             obj = obj.framesOverlap(Frame,index);%used only when
%             necessary
            obj = obj.pointsRFandCF();
            obj = obj.TempTT();
            obj = obj.meanTemperatureTool();
            obj = obj.maxTemperatureTool();
            obj = obj.maximumTemperature();
            obj = obj.maxTemperatureChip();
            obj = obj.calculateGradient();
            obj = extremePointsChip(obj);
            obj = obj.heatBalance(tuc,Vp,ap,widthTool);
            obj = obj.internalEnergyTool(widthTool);
            obj = obj.shearLine();
            obj = obj.calculatePecletNumber();
            obj = obj.forcesValues(Fp,Fq,ap,tuc);
            obj.TotalPowerBalance = 0.97*(obj.CuttingVelocity*(obj.CuttingForcePowerDirection*(1-A) + obj.CuttingForceParallelToolFace*A*obj.RatioR)/60);
            obj.HeatFluxThroughWorkpiece = obj.TotalPowerBalance - obj.HeatCarriedAwayByChip - obj.HeatFluxAwayFromToolTip;
        end
        
        function obj = framesOverlap(obj,Frame,index)%Used when it is necessary
            %to work with different emissivities
            cTT = obj.CoordinateToolTip;
            alpha = (90 - obj.ClearanceAngle)*pi/180;
            gamma = obj.RakeAngle*pi/180;
            p1 = cTT + 67*[-cos(gamma) sin(gamma)];
            p2 = cTT + 33*[-cos(alpha) sin(alpha)];
            c = [cTT(1) p1(1) p2(1)];
            r = [cTT(2) p1(2) p2(2)];
            biTool70 = roipoly(Frame(index).e70,c,r);
            aux = biTool70 == 1 & obj.biImageChip == 1;
            biTool70 = biTool70 - aux;
            biTool70andChip = biTool70 == 1 | obj.biImageChip == 1;
            biFrame85 = ones(size(Frame(index).e85)) - biTool70andChip;
            obj.frame = biTool70andChip.*Frame(index).e70 + biFrame85.*Frame(index).e85;
        end
        
        function obj = toolContour(obj)%Gets a binary image of the tool
            %Provide a binary image of the tool based on the coordinates of
            %the tool tip and on the angles of the faces (rake and
            %clearance faces)
            A = round(obj.CoordinateToolTip);
            m = size(obj.frame,1);
            xt = A(1);
            yt = A(2);
            y1 = round(yt + (xt - 1)*tan(obj.RakeAngle*pi/180));
            x2 = round(xt - (m - yt)*tan(pi/2 - (90 - obj.ClearanceAngle)*pi/180));
            c = [xt 0 0 x2];
            r = [yt y1 m m];
            B = roipoly(obj.frame,c,r);
            obj.biImageTool = B;
        end
        
        function obj = chipContour(obj)
            %method 1------------------------------------------------------
            %This method extracts the chip contour by taking the contour of
            %the line with 200 �C and then removes the tool contour(remaing
            %only the chip contour)
            c = obj.line200(1,:);
            r = obj.line200(2,:);
            B = roipoly(obj.frame,c,r);
            obj.biImageChip = B;
            B2 = obj.biImageTool ==1 & B == 1;
            B = B - B2;
            obj.biImageChip = B;
            imshow(B)
            K = B.*obj.frame;
            figure
            imagesc(K)
            %method 2------------------------------------------------------
            %This method is based on colors of the temperatures
            aux = round(obj.frame);
            aux(aux < 0) = 0;
            y = label2rgb(aux,'parula');
            image(y)
            lowerlimit = 175;
            if obj.UnCutChipThickness < 400*10^-6
                lowerlimit = 165;
            end
            B1 = y(:,:,2) > lowerlimit & y(:,:,3) < 210;
            B2 = obj.biImageTool ==1 & B1 == 1;
            B = B1 - B2;
            obj.biImageChip = B;
            figure
            imshow(B)
            K = B.*obj.frame;
            figure
            imagesc(K)
            figure
            imagesc(obj.frame)
        end
        
        function obj = maximumTemperature(obj)
            obj.MaximumTemperatureCuttingZone = max(max(obj.frame));
            [~,lin] = max(obj.frame);
            [~,col] = max(max(obj.frame));
            lin = lin(col);
            obj.xyMaxTemp = [col lin];%Position of the maximum temperature 
            %in the cutting zone
            imagesc(obj.frame)
            hold on
            plot(col,lin,'xr')
            hold off
        end
        
        function l = contactLength(obj)
            %Command to measure the length(in pixels) of the contact
            %between tool and chip
            imagesc(obj.frame)
            imdistline%Help to measure the amount of pixels on the contact length
            v = input('What is the value of the contact length for this frame [pixel]? ');
            close all
            l = 15*10^-6*v;
        end
        
        function obj = maxTemperatureTool(obj)
            %Gets the maximum temperature inside the tool contour
            C = obj.biImageTool;
            Frame = C.*obj.frame;
            T = max(max(Frame));
            obj.MaximumTemperatureTool = T;
        end
        
        function obj = maxTemperatureChip(obj)
            %Gets the maximum temperature inside the chip contour
            Frame = obj.biImageChip.*obj.frame;
            obj.MaximumTemperatureChip = max(max(Frame));
        end
        
        function obj = meanTemperatureTool(obj)
            %Get the points of valid pixels inside the tool and calculates
            %the average value for these points
            B = obj.biImageTool;
            Frame = B.*obj.frame;
            B = Frame > obj.validTemperature;
            Frame = B.*Frame;
            s = sum(sum(Frame));%Total sum of temperatures
            n = sum(sum(B));%Total number of valid pixels
            meanT = s/n;
            obj.MeanTemperatureTool = meanT;
        end       
        
        function obj = displayBinary(obj)%Plot the binary image of the frame
            imshow(obj.BW);
            hold on
            %Plot the points that determine the CF and RF and the tool tip
            %position
            plot(obj.coordRF(:,1),obj.coordRF(:,2),'bx')
            plot(obj.coordCF(:,1),obj.coordCF(:,2),'yx')
            plot(obj.CoordinateToolTip(1),obj.CoordinateToolTip(2),'xm')
            hold off
        end
        
        function obj = TempTT(obj)%Calculates the temperature at tool tip
            p1 = round(obj.CoordinateToolTip + 5*[-cos(obj.RakeAngle*pi/180) sin(obj.RakeAngle*pi/180)]);
            p2 = round(obj.CoordinateToolTip + 5*[-cos((90 - obj.ClearanceAngle)*pi/180) sin((90 - obj.ClearanceAngle)*pi/180)]);
            p3 = round(obj.CoordinateToolTip + 5*[-(cos(obj.RakeAngle*pi/180)+cos((90 - obj.ClearanceAngle)*pi/180)) (sin(obj.RakeAngle*pi/180)+sin((90 - obj.ClearanceAngle)*pi/180))]);
            T1 = obj.frame(p1(2),p1(1));
            T2 = obj.frame(p2(2),p2(1));
            T3 = obj.frame(p3(2),p3(1));
            TT = obj.frame(round(obj.CoordinateToolTip(2)),round(obj.CoordinateToolTip(1)));
            T = [T1 T2 T3 TT];
            obj.TemperatureToolTip = mean(T);
            %             imagesc(obj.frame)
            %             hold on
            %             plot([p1(1) p2(1) p3(1) obj.CoordinateToolTip(1)],[p1(2) p2(2) p3(2) obj.CoordinateToolTip(2)],'xm')
        end
        
        function obj = calculateCoordinates(obj)
            obj.BW = edge(obj.frame,'sobel');
            imshow(obj.BW)
            figure
            imagesc(f.*(ones(size(f))-obj.BW))
            %-------------------------Finding the clearance face-----------
            [H, THETA, RHO] = hough(obj.BW,'Theta',2:5);%Hough transformation
            %used to find angles between 2 - 5
            P  = houghpeaks(H, 10);
            obj.lines = houghlines(obj.BW, THETA, RHO, P, 'FillGap', 15,'MinLength',10);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
            l = length(obj.lines);
            obj.coordCF = [];
            for i=1:l%Loop to test the resultant lines from the command houghlines
                Theta=obj.lines(i).theta;
                t1 = obj.lines(i).point1;
                t2 = obj.lines(i).point2;
                rho = obj.lines(i).rho;
                imagesc(f)
                hold on
                plot([t1(1) t2(1)]',[t1(2) t2(2)]','m')
                hold off
                if rho < 204 && rho > 198
                    obj.coordCF = [t1;t2];
                    obj.ClearanceAngle = Theta;
                end
            end
            %-------------------------Finding the rake face----------------
            [H, THETA, RHO] = hough(obj.BW,'Theta',81:85);%Hough transformation
            %used to find angles between 5 - 9 (must be the complementary on this direction)
            P  = houghpeaks(H, 10);
            obj.lines = houghlines(obj.BW, THETA, RHO, P, 'FillGap', 15,'MinLength',10);%Here we can find the lines of cutting edge and afterwards find the coordinate of the tool tip
            l = length(obj.lines);
            obj.coordRF = [];
            for i=1:l
                Theta=obj.lines(i).theta;
                t1 = obj.lines(i).point1;
                t2 = obj.lines(i).point2;
                rho = obj.lines(i).rho;
                imagesc(f)
                hold on
                plot([t1(1) t2(1)]',[t1(2) t2(2)]','m')
                hold off
                if rho < 103 && rho > 98
                    obj.coordRF = [t1;t2];
                    obj.RakeAngle = 90 - Theta;
                end
            end
        end
        
        function obj = coordinateToolTip(obj)
            a = (obj.coordRF(1,2) - obj.coordRF(2,2))/(obj.coordRF(1,1) - obj.coordRF(2,1));%The slope of the rake face hardly will be Inf(Infinite) or NaN(Not-a-number),
            %because we took for this face a slope smaller than 45?
            b = obj.coordRF(1,2) - a*obj.coordRF(1,1);
            m = (obj.coordCF(1,2) - obj.coordCF(2,2))/(obj.coordCF(1,1) - obj.coordCF(2,1));%Slope of the cf, in some cases may be Inf(inclination of 90?, for example)
            h = @(x)(a*x + b);%line of the clearance face represented by f
            if m == Inf || m == -Inf%if the slope of the cf is 90? or -90?(Inf or -Inf)
                xi = obj.coordCF(1,1);%xi represents the coordinate x of the intersection(tool tip)
            else
                n = obj.coordCF(1,2) - m*obj.coordCF(1,1);
                xi = (n - b)/(a - m);
            end
            yi = h(xi);
            obj.CoordinateToolTip = [xi yi];
        end
        
        function obj = displayImageAndToolTip(obj)
            %Plots the fram image and the tool tip point
            figure
            imagesc(obj.frame);
            hold on
            plot(obj.CoordinateToolTip(1),obj.CoordinateToolTip(2),'xm')
            hold off
        end
        
        function obj = pointsRFandCF(obj)
            %This function gets points to help to plot temperature behavior
            %along the rake, clearance face and the bisectrix between these
            %faces
            alpha = (90 - obj.ClearanceAngle)*pi/180;
            gamma = obj.RakeAngle*pi/180;
            obj.pointRF = obj.CoordinateToolTip + 90*[-cos(gamma) sin(gamma)];
            obj.pointCF = obj.CoordinateToolTip + 90*[-cos(alpha) sin(alpha)];
            obj.pointM = obj.CoordinateToolTip + 40*[-2*cos(alpha)-cos(gamma) 2*sin(alpha)+sin(gamma)];
        end
        
        function vT = temperatureRFandCF(obj)
            pixelpitch = 15*10^-3;% mm/pixel
            extCF = obj.pointCF;% final point on the clearance face
            extRF = obj.pointRF;% final point on the rake face
            extM = obj.pointM;%final point of the bissectrix
            l1 = round(abs(obj.CoordinateToolTip(1)-extRF(1)));%length in pixels of rake line
            l2 = round(abs(obj.CoordinateToolTip(2)-extCF(2)));%length in pixels of clearance line
            l3 = max(round(abs(obj.CoordinateToolTip-extM)));
            vRFx = round(linspace(obj.CoordinateToolTip(1),extRF(1),l1));%coordinates x of the rake line
            vRFy = round(linspace(obj.CoordinateToolTip(2),extRF(2),l1));%coordinates y of the rake line
            vCFx = round(linspace(obj.CoordinateToolTip(1),extCF(1),l2));%coordinates x of the clearance line
            vCFy = round(linspace(obj.CoordinateToolTip(2),extCF(2),l2));%coordinates y of the clearance line
            vMx = round(linspace(obj.CoordinateToolTip(1),extM(1),l3));
            vMy = round(linspace(obj.CoordinateToolTip(2),extM(2),l3));
            T_RF = zeros(1,l1);%temperature for each pixel (each coordinate pair) - rake line
            T_CF = zeros(1,l2);%temperature for each pixel (each coordinate pair) - clearance line
            T_M = zeros(1,l3);
            for t=1:l1
                T_RF(t) = obj.frame(vRFy(t),vRFx(t));%Building the temperature vector - rake line
            end
            for t=1:l2
                T_CF(t) = obj.frame(vCFy(t),vCFx(t));%Building the temperature vector - clearance line
            end
            for t=1:l3
                T_M(t) = obj.frame(vMy(t),vMx(t));%Building the temperature vector - clearance line
            end
            d1 = zeros(1,l1);%distance for each pixel along the line
            d2 = zeros(1,l2);
            d3 = zeros(1,l3);
            for t=1:l1 - 1
                d1(t+1)=(((vRFx(t+1)-vRFx(1))^2)+((vRFy(t+1)-vRFy(1))^2))^(1/2);
            end
            for t=1:l2 - 1
                d2(t+1)=(((vCFx(t+1)-vCFx(1))^2)+((vCFy(t+1)-vCFy(1))^2))^(1/2);
            end
            for t=1:l3 - 1
                d3(t+1)=(((vMx(t+1)-vMx(1))^2)+((vMy(t+1)-vMy(1))^2))^(1/2);
            end
            d1 = d1*pixelpitch;
            d2 = d2*pixelpitch;
            d3 = d3*pixelpitch;
            figure
            hold on
            plot(d1,T_RF)
            plot(d2,T_CF)
            plot(d3,T_M)
            xlabel('Distance from the tool tip (mm)')
            ylabel('Temperature (�C)')
            legend('Rake face','Clearance face','Middle vector')
            %saveas(fig,file,'jpeg')
            hold off
            figure
            imagesc(obj.frame)
            colormap jet
            hold on
            plot(vRFx,vRFy,'k','LineWidth',1)
            plot(vCFx,vCFy,'k','LineWidth',1)
            plot(vMx,vMy,'k','LineWidth',1)
            hold off
            m = min([l1 l2 l3]);
            vT = [d1(1:m)' T_RF(1:m)' d2(1:m)' T_CF(1:m)' d3(1:m)' T_M(1:m)'];
            %Not only plot the graph but also export the points of
            %Temperature x distance
        end
        
        function obj = extremePointsChip(obj)
            [y,x] = find(obj.lineChip);%gets a range of the coordinates of this line
            %The first point is on the highest position and the last is on
            %the lowest(which are the extremity of the line)
            obj.extPtosLineChip = [x(1) y(1);x(end) y(end)];
        end
        
        function obj = displayIsotherms(obj)
            tRF = obj.RakeAngle*pi/180;
            tCF = (90 - obj.ClearanceAngle)*pi/180;
            vRF = [-cos(tRF) sin(tRF)];
            vCF = [-cos(tCF) sin(tCF)];
            %p1 RF direction
            t = (obj.CoordinateToolTip(1) - 1)/vRF(1);
            p1 = obj.CoordinateToolTip - t*vRF;
            %p2 CF direction
            t = (256 - obj.CoordinateToolTip(2))/vCF(2);
            p2 = obj.CoordinateToolTip + t*vCF;
            %auxiliar to plot
            auxX = [p1(1) obj.CoordinateToolTip(1) p2(1)]';
            auxY = [p1(2) obj.CoordinateToolTip(2) p2(2)]';
            Tmax = max(max(obj.biImageTool.*obj.frame));
            Tv = obj.validTemperature;
            v = round(Tv:40:Tmax);
            %Display only isotherms and tool contour-----------------------
            figure
            [C,h] = contour(obj.frame,v);
            hold on
            clabel(C,h,'FontSize',10)
            h.LineWidth = 1.5;
            daspect([1,1,1])
            plot(auxX,auxY,'k')
            text(150,150,'Tool')
            text(100,75,'Chip')
            xlabel('pixel');
            ylabel('pixel');
            title('Isotherms');
            cb = colorbar('vert');
            zlab = get(cb,'ylabel');
            set(zlab,'String','Temperature (�C)');
            axis([50 250 25 200])
            %saveas(fig,t,'jpeg')
            axis off
            colormap jet
            hold off
            %Display tool and isotherms-----------------------------------
            lc = obj.extPtosLineChip;
            figure
            imagesc(obj.frame)
            colormap jet
            hold on
            plot(auxX,auxY,'k')
            plot(lc(:,1),lc(:,2),'k--','LineWidth',1)%line in the chip on
            %the end of contact
            [C,h] = contour(obj.frame,v);
            h.LineColor = [0.247 0.247 0.247];%Set only one color for the
            %contour lines
            clabel(C,h,'manual','FontSize',10);
            x = obj.CoordinateToolTip(1);
            y = obj.CoordinateToolTip(2);
            axis([x-180 x+15 y-60 y+110])%resize the image
            daspect([1,1,1])
            axis off
            hold off
        end
        
        function obj = calculateGradient(obj)
            pp = 15*10^-6;%pixelpitch in m
            tx = zeros(size(obj.frame));
            ty = zeros(size(obj.frame));
            k = 0;
            for j = 1:5%loop to calculate the gradients of temperature
                [auxx,auxy]=gradaux_v2(obj.frame,j);
                tx = tx + auxx;
                ty = ty + auxy;
                k = k + 1;
            end
            obj.Tx = tx/(k*pp);
            obj.Ty = ty/(k*pp);
        end
        
        function obj = displayGradient(obj)
            %Display the gradient arrows for tool tip are, part of the rake
            %face and clearance face
            auxx = [obj.pointCF(1) obj.CoordinateToolTip(1) obj.pointRF(1)];
            auxy = [obj.pointCF(2) obj.CoordinateToolTip(2) obj.pointRF(2)];
            k = obj.heatConductivity;
            qx = -k*obj.Tx;
            qy = -k*obj.Ty;
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 10;
            xmax = obj.CoordinateToolTip(1) + 5;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 10;
            axis([xmin xmax ymin ymax])
            title('Tool Tip')
            daspect([1,1,1])
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 30;
            xmax = obj.CoordinateToolTip(1) - 10;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 15;
            axis([xmin xmax ymin ymax])
            title('Rake Face')
            daspect([1,1,1])
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            xmin = obj.CoordinateToolTip(1) - 10;
            xmax = obj.CoordinateToolTip(1) + 10;
            ymin = obj.CoordinateToolTip(2) + 10;
            ymax = obj.CoordinateToolTip(2) + 20;
            axis([xmin xmax ymin ymax])
            title('Clearance Face')
            daspect([1,1,1])
        end
        
        function obj = displayGradientContour(obj)
            auxx = [obj.pointCF(1) obj.CoordinateToolTip(1) obj.pointRF(1)];
            auxy = [obj.pointCF(2) obj.CoordinateToolTip(2) obj.pointRF(2)];
            k = obj.heatConductivity;
            qx = -k*obj.Tx;
            qy = -k*obj.Ty;
            figure
            quiver(qx,qy)
            hold on
            plot(auxx,auxy,'k')
            contour(obj.frame,10)
            xmin = obj.CoordinateToolTip(1) - 20;
            xmax = obj.CoordinateToolTip(1) + 5;
            ymin = obj.CoordinateToolTip(2) - 5;
            ymax = obj.CoordinateToolTip(2) + 20;
            axis([xmin xmax ymin ymax])
            daspect([1,1,1])
        end
        
        function obj = findLineChip(obj)%Gets the line perpendicular to the
            %chip movement on the end of the contact tool - chip
            [m,n] = size(obj.frame);
            o = obj.RakeAngle*pi/180;
            l = obj.ContactLength/(15*10^-6);
            c = obj.CoordinateToolTip + l*[-cos(o) sin(o)];
            xm = c(1);
            ym = c(2);
            x1 = xm - tan(o)*(ym - 1);
            x2 = x1 + tan(o)*(m - 1);
            vx = round(linspace(x1,x2,m));
            vy = linspace(1,m,m);
            B1 = zeros(m,n);
            for i = 1:m
                B1(vy(i),vx(i)) = 1;
            end
            B2 = B1 == 1 & obj.biImageChip == 1;
            obj.lineChip = B2;
            imshow(B2)
            figure
            imshow(obj.biImageChip)
            aux = (obj.biImageChip - B2).*obj.frame;
            figure
            imagesc(aux)
        end
        
        function obj = findLineTool(obj)%Function that gets the coordinates of
            %the contour lines for each temperature level
            [m,n] = size(obj.frame);
            Tmax = max(max(obj.biImageTool.*obj.frame));
            Tv = obj.validTemperature;
            v = round(Tv:40:Tmax);%the step between the temperatures is 40
            if length(v) == 1%Conditional if the min and max temperatures are equal
                v = round([Tv Tmax]);
            end
            [C,~] = contour(obj.frame,v);%Getting the contour lines for C
            close
            l = length(v);
            B = zeros(m,n,l);
            C = round(C);
            for k = 1:l
                [~,J] = find(C == v(k));
                [~,p] = max(C(2,J));
                J = J(p);
                for z = J+1:J+C(2,J)
                    B(C(2,z),C(1,z),k) = 1;
                end
                if k == 1%Take the last line for the future heat flow calculation
                    obj.line200 = C(:,J+1:J+C(2,J));
                end
                B(:,:,k) = B(:,:,k).*obj.biImageTool;
            end
            obj.lineTool = B;
            for i = 1:l %Loop to check the binary image of the lines
                imshow(B(:,:,i))
            end
        end
        
        function obj = heatBalance(obj,tuc,Vc,ap,w)
            pp = 15*10^-6;% pixel pitch [m/pixel]
            %-------------------------------------------------------------------------
            %First part - Heat carried away by the chip
            k = obj.heatConductivity;% Heat conductivity tool
            cp = obj.heatCapacityWP;
            M = obj.lineChip.*obj.frame;
            MH = polyval(cp,M);
            MH(MH == cp(end)) = 0;%Remove the values for tempertures equal to zero
            Ht = MH.*(obj.frame-22);%J/kg - 22 is the temperature of the environment
            Ht = sum(sum(Ht));
            n = sum(sum(obj.lineChip));
            Hc = Ht/n; %mean entalpy on the line chip
            p = 7697; %density of workpiece - kg/m^3
            Qc = Hc*Vc*tuc*p;%Vc*tuc is equals to Vchip*tchip
            obj.HeatCarriedAwayByChip = Qc*ap;
            %-------------------------------------------------------------------------
            %Second part - Heat carried away by the tool
            dT = ((obj.Tx).^2 + (obj.Ty).^2).^(1/2);%Resultant gradient of temperature
            Q = zeros(size(obj.lineTool,3),1);
            for i = 1:size(obj.lineTool,3)
                L = obj.lineTool(:,:,i);
                Q(i) = sum(sum(L.*dT))*pp*k*w;
            end
            obj.heatAccumulatedPerLine = Q;
            Qm = mean(Q(1:2));%mean of the last two lines
            obj.HeatFluxAwayFromToolTip = Qm;
        end
        
        function n = exceedingPoints(obj,Temperature)
            B = obj.frame.*obj.biImageTool > Temperature;
            n = sum(sum(B));
%             K = sprintf('Number of points exceeding the defined temperature is: %3d \n',n);
%             disp(K)
        end
        
        function obj = internalEnergyTool(obj,w)%Calculates the internal energy
            %in the tool, but only for valid pixels
            pp = 15*10^-6;%pixel pitch in cm
            cp = obj.heatCapacityTool;%J/(K*m^3)
            %Heat capacity for tool
            Te = 22;%Environment
            B = obj.frame.*obj.biImageTool > obj.validTemperature;
            B1 = obj.frame.*B;
            B2 = polyval(cp,B1);%Heat capacity for each pixel (J/kgK)
            B2(B2 == cp(4)) = 0;
            H = B2.*(obj.frame - Te)*(pp^2);%Heat Amount for each pixel(J/m)
            Ha = sum(sum(H));%Mean value for the entire tool
            obj.InternalEnergyTool = Ha*w;%w is the width of the workpiece
        end
        
        function obj = shearLine(obj)%This function tries to identify the
            %of the primary shear zone (to get the shear angle)
            B = obj.biImageChip;
            v1 = sum(B);
            v1(v1 == 0) = [];
            l1 = length(v1);
            C = imcrop(B,[20 20 l1 100]);%Cutting the image to a new size
            %to focus on the shear region
            [m,n] = size(C);
            pto = zeros(1000,2);
            count = 1;
            %In the following loop we create a binary image of the chip
            %contour without its boundaries which is necessary to identify the shear line
            for i = 2:m-1
                for j = 2:n-1
                    if C(i,j+1) == 1 && C(i,j-1) == 1 && C(i+1,j) == 1 && C(i-1,j) == 1
                        pto(count,:) = [i j];
                        count = count + 1;
                    end
                end
            end
            %On the following loop the empty cells are removed
            for i =1000:-1:1
                if isequal(pto(i,:),[0 0]) == 1
                    pto(i,:) = [];
                end
            end
            l = size(pto,1);
            %In this loop only the boundaries of the chip is left
            for i = 1:l
                C(pto(i,1),pto(i,2)) = 0;
            end
            
            [H, THETA, RHO] = hough(C,'Theta',-40:-30);%Hough transformation
            %from 30 until 40 is the range of possible angles for the shear
            %angle
            P  = houghpeaks(H, 5);
            lin = houghlines(C, THETA, RHO, P, 'FillGap', 15,'MinLength',10);
            l=length(lin);
            %p1 and p2 will be the extreme points of the shear line
            p1 = [];
            p2 = [];
            for i=1:l
                Theta=lin(i).theta;
                t1 = lin(i).point1;
                t2 = lin(i).point2;
                y = abs(t1(2)-t2(2));
%                 imshow(C)
%                 hold on
%                 plot([t1(1) t2(1)]',[t1(2) t2(2)]','xm')
%                 hold off
                if isempty(p1) && isempty(p2) && abs(Theta + 34) < 5
                    p1 = t1 + [19 19];%Plus 19 due to the previous cutting image process
                    p2 = t2 + [19 19];
                    ym = y;
                    obj.ShearAngle = abs(Theta);
                end
                if abs(Theta + 34) < 5 && y > ym
                    p1 = t1 + [19 19];
                    p2 = t2 + [19 19];
                    obj.ShearAngle = abs(Theta);
                end
            end
            if isempty(obj.ShearAngle)%It was not possible to find a good line
                obj.ShearAngle = 30;%Default condition
            end
        end
        
        function obj = forcesValues(obj,Fp,Fq,ap,tuc)
            phi = obj.ShearAngle*pi/180;%shear angle
            gamma = obj.RakeAngle*pi/180;%Rake angle
            Fs = Fp*cos(phi) - Fq*sin(phi);%Cutting force component parallel to shear plane
            Ns = Fq*cos(phi) + Fp*sin(phi);%Cutting force component perpendicular to shear plane
            Fc = Fp*sin(gamma) + Fq*cos(gamma);%Cutting force component parallel to tool face
            Nc = Fp*cos(gamma) - Fq*sin(gamma);%Cutting force component perpendicular to tool face
            mu = Fc/Nc; % coefficient of friction
            As = ap*tuc/sin(phi);%Area shear plane
            tau = Fs/As;%shear stress
            sigma = Ns/As;%Normal stress
            r = sin(phi)/cos(phi - gamma); %ratio r = t/tc = lc/l
            ss = cos(gamma)/(sin(phi)*cos(phi-gamma));%shear strain
            us = tau*ss;%shear energy per volume
            uf = Fc*r/(tuc*ap);%friction energy per volume
            beta = atan(Fc/Nc);%friction angle on tool face
            obj.CuttingForceParallelToolFace = Fc;
            obj.CuttingForcePowerDirection = Fp;
            obj.CuttingForceUncutChipThicknessDirection = Fq;
            obj.CuttingForceParallelShearPlane = Fs;
            obj.CuttingForcePerpendicularShearPlane = Ns;
            obj.CuttingForcePerpendicularToolFace = Nc;
            obj.CoefficientFriction = mu;
            obj.ShearStress = tau;
            obj.NormalStress = sigma;
            obj.RatioR = r;
            obj.ShearEnergyVolume = us;
            obj.FrictionEnergyVolume = uf;
            obj.FrictionAngle = beta*180/pi;
        end
        
        function obj = calculatePecletNumber(obj)
            cp = polyval(obj.heatCapacityWP,obj.MaximumTemperatureCuttingZone);
            %heat capacity for the maximum temperature(J/kgK)
            k = 39.1;%heat conductivity of the workpiece (W/mK) - mean value
            d = 7.697*10^3;%density of the tool (kg/m^3)
            obj.PecletNumber = ((obj.CuttingVelocity/60)*obj.UnCutChipThickness)/(k/(cp*d)); 
        end
    end
end