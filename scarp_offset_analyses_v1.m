
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      Christine Regalla %
%                                                      last modified     %
%                                                      3/23/16           %
%                                                                        %     
%                                                                        %
% This script will calcualte vertical separation, heave, throw, and      %
% fault slip from a topographic profile across a fault scarp. This       %
% code reads a tab demited text file of topographic profile data, which  %
% must be formatted in two columns of x (distance) and z (elevation)     %
% data points. Once text profile is entered, the profile can be saved as % 
% a .mat file and reloaded in subsequent runs. Midpoint and regressions  %
% through lower and upper surfaces can be chosen graphically or entered  %
% numerically. Requires Matlab Statistics Toolbox.                       %
%                                                                        %
% Inputs:                                                                %
%   filename.txt = tab delimited text file (without formatting)          %
%           containing x and z data for the topographic scarp profile    %
%                                                                        %
% Outputs:                                                               %
%    VS = vertical separation, calculated as the difference between the  %
%          elevation of the projection of upper and lower surfaces at    %
%           the scarp midpoint                                           %
                                                    %
%    ru = slope and intercept of regression through lower surface        %
%    rl = slope and intercept of regression through upper surface        %
%    ru_l = intersection point between lower fan surface and scarp face  %
%    rl_i = intersection point between upper fan surface and scarp face  %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;  %clear all existing variables; close all existing plots

         %___________________________________________________%

         % SET SCRIPT VARIABLES %
              
             FAULTANGLE = 60;    %fault angle (degrees)  positive (+) faults dips have a slope that is 
                                    % positive in x,y, space; that is the fault plane slopes down to the left
                                    % toward the negative x axis                      
             FAULTUNCERT = 20;      % uncertainty (degree) in fault anlge
             dx = 0.25;          % x-interval used to discretize profile
                 
             nruns = 100;        % # of Monte Carlo Runs
                          
             alpha = .05;       %confidence interval for regression uncertainty (e.g. .05 for 95%...2sigma, or .01 for 99%....3sigma)
 
         %___________________________________________________%
            
                %Convert FAULTANLGE to dimensionless slope 
                 mf = tan(FAULTANGLE*pi/180);
                 
                 
     % READ IN INPUT TEXT FILE FILENAME.txt     
     qq = (input('Import a) raw profile data (.txt)? or b) processed data file (.mat)? (a or b):     ', 's'));
     filename = (input('Enter the profile filename, with extension:     ', 's'));
     if qq == 'a'
        p=dlmread(filename);  
        plot (p(:,1), p(:,2), '.b');

    
    %REVERSE ORIENTATION OF PROFILE IF RIGHT-FACING
        o = (input('Does the scarp face a)left or b)right? (a or b):     ', 's')); 
        if o=='b'
            p2= [[p(2:end,1)-p(1:end-1,1);0], p(:,2)];
             for i=1:size(p,1)
                p3(i,:) = p2(end-i+1, :);        
             end
            p2(:,:)=[cumsum(p3(:,1)), p3(:,2)];
            p=p2;
        end
        clf

    
    %DISCRETIZE PROFILE%
    x=[ceil(p(1,1)):dx:floor(p(end,1))]';   
    z = interp1(p(:,1), p(:,2), x); %interpolate initial profile over dx
    hold on; 
    h(1) = plot (p(:,1), p(:,2), '.b'); h(2) = plot (x,z, '-k');
    xlabel 'Distance Along Profile'; ylabel 'Height';
    legend (h, 'Profile points', 'Interpolated Profile', 'Location','SouthEast');
    
    
    %SELECT MIDPOINT AND REGRESS BOUNDING SURFACES
    q = (input('Do you wish to select midpoint and upper, lower surfaces a) graphically or b) by entering numerical values? (a or b):     ', 's'));
    h=1; hu=1; hl= 1;
        %SELECT MIDPOINT%
        a='n';
        while a=='n' 
            if q=='a'
                disp(' ')
                disp('Select the midpoint of the profile.     ')
                     
        %{                
        % workaround text to make the ginput cursur to appear - a bug in old
        % Matlab version
            %pcolor (64*rand(100));

        legend('location','westoutside')

        legend('hide')
        %}
                
                midpt = ginput(1); 
            elseif q=='b'
                midpt = input('Enter midpoint:  [x,z]    ');
            end
            d = sqrt((x-midpt(1)).^2 + (z-midpt(2)).^2);
            midpt = [x(d==min(d)), z(d==min(d))];
            temppoint = plot(midpt(1), midpt(2), 'or');
            disp(' '); a=input('Would you like to keep the current midpt? (y/n)    ', 's');
            if a=='n'
                delete(temppoint)
            end
        end
        
        
        %SELECT LOWER AND UPPER REGRESSION SURFACES%
        a='n';
        while a=='n' 
            if q=='a'
                disp(' ')
                disp('Select range of points to use for lower surface regression, from left to right.     ')
                rl_range = ginput(2); 
            elseif q=='b'
                rl_range = input('Enter bounding points on lower surface:  [left_x, left_z; right_x, right_z]    ');   
            end
            %delete (hl);
            d = sqrt((x-rl_range(1,1)).^2 + (z-rl_range(1,2)).^2);
            ind = find (d==min(d));
            d = sqrt((x-rl_range(2,1)).^2 + (z-rl_range(2,2)).^2);
            ind(2) = find (d==min(d));   
            rl = polyfit(x(ind(1):ind(2)), z(ind(1):ind(2)), 1);
   
            hl = refline(rl);
            disp(' '); a=input('Would you like to keep the regression? (y/n)    ', 's');
            if a=='n'
                delete(hl)
            end
            
        end
            rl_range = [x(ind(1):ind(2)), z(ind(1):ind(2))];
        
        a='n';
        while a=='n'  
            if q=='a'
                disp(' ')
                disp('Select range of points to use for upper surface regression, from left to right.     ')
                ru_range = ginput(2); 
            elseif q=='b'
                ru_range = input('Enter bounding points on upper surface:  [left_x, left_z; right_x, right_z]    ');
            end
            d = sqrt((x-ru_range(1,1)).^2 + (z-ru_range(1,2)).^2);
            ind = find (d==min(d));
            d = sqrt((x-ru_range(2,1)).^2 + (z-ru_range(2,2)).^2);
            ind(2) = find (d==min(d));
            ru = polyfit(x(ind(1):ind(2)), z(ind(1):ind(2)), 1);
         
            hu = refline(ru);
            disp(' '); a=input('Would you like to keep the regression? (y/n)    ', 's');
            if a=='n'
                delete(hu)
            end
        end
            ru_range = [x(ind(1):ind(2)), z(ind(1):ind(2))];
      
            
        %SELECT RANGE OF COORDINATES WHERE FAULT PLANE MAY INTERSECT SCARP FACE%              
        a='n';
        while a=='n'  
            if q=='a'
                disp(' ')
                disp('Select range of points that bound the possible intersection of the fault plane with the scarp face, from left to right.     ')
                rf_range = ginput(2);
            elseif q=='b'
                rf_range = input('Enter bounding points on scarpface:  [left_x, left_z; right_x, right_z]    ');
            end
                d = sqrt((x-rf_range(1,1)).^2 + (z-rf_range(1,2)).^2);
                ind = find (d==min(d));
                d = sqrt((x-rf_range(2,1)).^2 + (z-rf_range(2,2)).^2);
                ind(2) = find (d==min(d));
                rfs = polyfit(x(ind(1):ind(2)), z(ind(1):ind(2)), 1);
                
                rf_range = [x(ind(1):ind(2)), z(ind(1):ind(2))]; 
                
                
                hf = plot(rf_range(:,1), rf_range(:,2), '-r', 'LineWidth', 2);
                disp(' '); a=input('Would you like to keep this range? (y/n)    ', 's');
            if a=='n'
                delete(hf)
                clear rf_range
            end
        end
                 
               
            
        %CENTER PROFILE    
     
        
        x=x-midpt(1); z = z-midpt(2);
        
        rl_range(:,1) = rl_range(:,1) - midpt(1); rl_range(:,2) = rl_range(:,2) - midpt(2); 
        ru_range(:,1) = ru_range(:,1) - midpt(1); ru_range(:,2) = ru_range(:,2) - midpt(2); 
        rf_range(:,1) = rf_range(:,1) - midpt(1); rf_range(:,2) = rf_range(:,2) - midpt(2); 

            % Regress:
            % [b,bint] = regress(y,x); 
            % b = [y-intercept, slope], best value
            % bint = [min y-int, max y-int; min slope, max slope];
        
        
            % Recalc y=mx+b for Upper surface  
                    x1 = ru_range(:,1);  z1 = ru_range(:,2);
                    X1 = [ones(size(x1)) x1];
                    [temp temp_int] = regress(z1,X1);
                    ru(1)=temp(2); ru(2) = temp(1);
                    ru_int(1,:) = temp_int(2,:); ru_int(2,:) = temp_int(1,:);
        
            % Recalc y=mx+b for Lower surface
                    x1 = rl_range(:,1);  z1 = rl_range(:,2);
                    X1 = [ones(size(x1)) x1];
                    [temp temp_int] = regress(z1,X1);
                    rl(1)=temp(2); rl(2) = temp(1);
                    rl_int(1,:) = temp_int(2,:); rl_int(2,:) = temp_int(1,:); 
        
            % Recalc y=mx+b for Fault scarp face:
                    x1 = rf_range(:,1);  z1 = rf_range(:,2);
                    X1 = [ones(size(x1)) x1];
                    [temp temp_int] = regress(z1,X1);
                    rfs(1)=temp(2); rfs(2) = temp(1);
                    %rfs_int(1,:) = temp_int(2,:); rfs_int(2,:) = temp_int(1,:);

        
     % PLOT SCARP
        clear h
        p(:,1) = p(:,1) - midpt(1); 
        p(:,2) = p(:,2) - midpt(2); 
        clf; hold on; 
        h(1) = plot (x,z, '-k', 'LineWidth', 2); h(2) = refline(ru); refline(rl); 
            xlabel ('Distance (m)')
            ylabel ('Elevtion (m)')
         
         %plot scarp face regression 
         %hf = plot(rf_range(:,1), rf_range(:,2), '-r', 'LineWidth', 2);

       %
         % Plot confidence intervals on upper surface:
             % http://www.real-statistics.com/regression/confidence-and-prediction-intervals/
             zhat = ru(2)+ru(1)*x;   %calcualte z values for regression line
             t = tinv(1-alpha/2,(length(x)-2)); %t value of distribution
             z_avg = mean(z); x_avg = mean(x);
             SX = sum((x-x_avg).^2); SY = sum((z-z_avg).^2);
             SXSY = sum((x-x_avg).*(z-z_avg)).^2;
             syx = sqrt((1/(length(x)-2))*(SY-(SXSY/SX)));  %std error in estimate ... Syx ... primary equation to solve
             CI = t*syx*sqrt(1/(length(x))+(x(:)-x_avg).^2/SX); % Calculate Confidence interval

             %calculate upper and lower confidence curves (linear)
             z_plus_CI = zhat + CI;
             z_minus_CI = zhat - CI;

             %plot confidence intervals
              plot(x,z_plus_CI,'r-.');
              plot(x,z_minus_CI,'r-.');
              
         % Plot confidence intervals on upper surface:
            % http://www.real-statistics.com/regression/confidence-and-prediction-intervals/
             
             zhat = rl(2)+rl(1)*x;   %calcualte z values for regression line
             t = tinv(1-alpha/2,(length(x)-2)); %t value of distribution
             z_avg = mean(z); x_avg = mean(x);
             SX = sum((x-x_avg).^2); SY = sum((z-z_avg).^2);
             SXSY = sum((x-x_avg).*(z-z_avg)).^2;
             syx = sqrt((1/(length(x)-2))*(SY-(SXSY/SX)));  %std error in estimate ... Syx ... primary equation to solve
             CI = t*syx*sqrt(1/(length(x))+(x(:)-x_avg).^2/SX); % Calculate Confidence interval

             %calculate upper and lower confidence curves (linear)
             z_plus_CI = zhat + CI;
             z_minus_CI = zhat - CI;


             %plot confidence intervals
              h(3) = plot(x,z_plus_CI,'r-.');
              plot(x,z_minus_CI,'r-.');
           
       %}
       
  %%% save %%%     
    disp(' ')    
    q = (input('Save processed profile as .mat file? (y/n)       ', 's'));
    if q=='y'
        name = [input('Enter file name:   ', 's'), '.mat'];
        save (name)
        %save (name, 'midpt', 'p', 'x', 'z', 'rl_range', 'rl', 'rl_int', 'ru_range', 'ru', 'ru_int', 'rfs', 'rf_range')
    end

    
% LOAD MAT FILE
% Load profile from Matlab variables if you are not laoding a text file

%%%%% PLOT CI %%%%%%

 elseif qq=='b'
     load (filename)
     clf; hold on; h(1) = plot (x,z, '-k', 'LineWidth', 2); h(2) = refline(ru); refline(rl);
             %plot confidence intervals
              h(3) = plot(x,z_plus_CI,'r-.');
              plot(x,z_minus_CI,'r-.');
     xlabel ('Distance (m)'); ylabel ('Elevation (m)')
       
 end

 
    
%%  MONTE CARO SIMULATION TO CALCULATE VERTICAL SEPARATION, SLIP, HEAVE, THROW 
     
    % Equation for upper surface is:  y = mu*x + bu
    % Equation for lower surface is:  y = ml*x + bu
    % Equation for fault plane is:  y = mf*x + bf
        % Equation for the plane regressed through the scarp face is:   
                                   % y = rf(1)*x + rf(2)
    
    % where:
        % mu = slope of upper surface;  ru(1)
        % bu = y-int of upper surface; ru(2)
        % ml = slope of lower surface;  rl(1)
        % bl = y-int of lower surface; rl(2)
        % mf = slope of fault plane (= tan(FAULTANGLE))
        % bf = y-int of fault plane (right now it is 0)
        
        
    % therefore:  x = ((bu-bf)/(mf-mu));
    % and:        y = mf*((bu-bf)/(mf-mu)) + bf;

  
   
 for i=1:nruns 
     
    % 1) Randomly generate fault angles
         F = (FAULTANGLE-FAULTUNCERT) + rand(1)* (2*FAULTUNCERT);
         FAULTANGLES(i) = F;
            mf = tan(F.*pi/180);
    % 2) Randomly generate a y-intercept for fault plane (given input range of coordiantes where fault can intersect the fault scarp)
        % a) select point on scarp face where fault can intersect (choose an x value, fsx)
            fsx = rf_range(1,1) + rand(1)*(rf_range(end,1)-rf_range(1,1));
        % b) solve for the y -coordinate on the scarp face regression, using rf = (mf, bf);
            fsy = rfs(1)*fsx+rfs(2);
                %plot (fsx, fsy, '.r', 'MarkerSize', 30);
        % c) use point slope formula to calculate y-int for the fault plane, yield bf
                % y - yf = mf*(x - xf)
                % y = mf*x + (-mf*fsx + fsy)
                % so b = fsy - mf*fsx
            bf = fsy - mf*fsx;
            
             
     % plot "fault"
            dy = 10;  %'depth" to which fault will be plotted
            h(4) = plot([fsx, fsx-((dy)/(mf))], [fsy, fsy-dy], '--k');
  
    % 3) Randomly generate a regression line for upper surface 
  
            %R = normrnd(mu,sigma)
            RU(1) = normrnd(ru(1), ((ru_int(1,2) - ru_int(1,1))/2));
            RU(2) = normrnd(ru(2), ((ru_int(2,2) - ru_int(2,1))/2));
            RUS(i,:) = RU;
      
           
     % 4) Randomly generate a regression line for lower surface
   
            RL(1) = normrnd(rl(1), ((rl_int(1,2) - rl_int(1,1))/2));
            RL(2) = normrnd(rl(2), ((rl_int(2,2) - rl_int(2,1))/2));
            RLS(i,:) = RL;
        
         
       
   % CALCULATE INTERSECTION BETWEEN FAULT PLANE AND UPPER / LOWER SURFACES         
     % P1 = intersection (x,y) of upper surface and fault plane
         P1 = [(RU(2)-bf)/(mf-RU(1)); ((RU(2)-bf)/(mf-RU(1)))*(mf)+bf];
     % P2 = intersection (x,y) of lower surface and fault plane
         P2 = [(RL(2)-bf)/(mf-RL(1)); ((RL(2)-bf)/(mf-RL(1)))*(mf)+bf];
         
  
    % Where:     
        % RU = a vector that contains the slope and intercept of the upper
                % surface regression line
        % RL = a vector that contains the slope and intercept of the lower
                % surface regression line
   
     % CALCUALTE THE INTERSECTION BETWEEN THE REGRESSION THROUGH THE FAULT
     % SCARP FACE AND THE UPPER / LOWER SURFACES
            % PA = intersection (x,y) of upper surface and scarp face
                PA = [(RU(2)-rfs(2))/((rfs(1))-RU(1)); (RU(2)-rfs(2))/(rfs(1)-RU(1))*rfs(1)+rfs(2)];
            % PB = intersection (x,y) of lower surface and scarp face
                PB = [(RL(2)-rfs(2))/((rfs(1))-RL(1)); (RL(2)-rfs(2))/(rfs(1)-RL(1))*rfs(1)+rfs(2)];

           % plot regression through scarp face
                %plot ([PA(1),PB(1)], [PA(2), PB(2)], '-b')
      
      

% CALCULATE VERTICAL SEPARATION                
% VS is therefore the difference in z values between where a veritical line
% through the intersection point of the fault and the scarp face (fsx, fsy), 
% intersects the upper and lower surfaces.  this is the difference between
% y(fsx) for the upper and lower surfaces
    
    VS(i) = (RU(1)*fsx+RU(2))-(RL(1)*fsx+RL(2)); 
        
    % Plot line along which vertical separation is measured
        h(5) = plot([fsx,fsx], [(RU(1)*fsx+RU(2)), (RL(1)*fsx+RL(2))], '-r');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULATE DIP SLIP ALONG FAULT % 
    fault_slip(i) = sqrt((P2(1,:)-P1(1,:)).^2 + (P2(2, :)-P1(2, :)).^2);  
        %plot slip vector
        % plot([P2(1), P1(1)],[P2(2), P1(2)], '-g' );


%
% CALCULATE HEAVE AND THROW
    %Heave = the difference between the x vlaues of the points of
    %intersection between the fault plane and the upper/ lower surfaces
        Heave(i) = abs(P1(1)-P2(1)); 
 
    %Throw = the difference between the y vlaues of the points of
    %intersection between the fault plane and the upper/ lower surfaces
        Throw(i) = abs(P1(2))+abs(P2(2)); 
%}
      
 end
  ave_vert_sep = mean(VS) 
    stdtev_vert_sep = std(VS)
 ave_slip = mean(fault_slip) 
    stdtev_slip = std(fault_slip)
 ave_heave = mean(Heave)  
    stdtev_heave = std(Heave)
 ave_throw = mean(Throw)    
    stdtev_throw = std(Throw) 
 
 plotlegend = legend([h], 'Topo Profile', 'Surface Regressions', '95% Regression CI',...
     'Fault Plane', 'Vertical Separation', 'Location', 'Northeast');
        


%%

%Save variables
disp(' ')    
q = (input('Save variables as .mat file, and plot as .fig? (y/n)       ', 's'));
if q=='y'
       corename = input('Enter file name:   ', 's');;
    name = [corename, '.mat'];
    save (name)
    save (name, 'midpt', 'p', 'x', 'z', 'rl_range', 'rl', 'ru_range', 'ru', 'rfs', 'rf_range',...
        'FAULTANGLES', 'VS','ave_vert_sep','stdtev_vert_sep', 'ave_slip', 'stdtev_slip', 'fault_slip',...
        'Heave', 'ave_heave', 'stdtev_heave', 'Throw', 'ave_throw', 'stdtev_throw');
    title(corename); saveas(gcf,corename);
end
    




    
