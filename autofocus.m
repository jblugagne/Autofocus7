function autofocus()
% This function uses the settings in guihs.autofocus.settings. You can swap 
% them for multiple positions


    global mmc guihs

    if isfield(guihs,'autofocus')
        dump = guihs.autofocus;
        if isfield(dump,'settings')
            AS = guihs.autofocus.settings;
        else
            AS = struct();
        end
    else
        AS = struct();
    end
    
    % Has setup been done already?
    if  ~(isfield(AS,'initialized') && AS.initialized)
        % If it has not been done, intialize:
        AS = init_autofocus();
    else
        % If it has, acquire sharpness curve: 
        initZ = mmc.getPosition('PIZStage');
        pause(.1) % What's this?
        
        % Calculate the range:
        nbsteps = round(.5*AS.CoarseR/AS.CoarseSS);
        switch AS.method
            case 'linear'
                Range = cumsum(repmat(AS.CoarseSS,1,nbsteps));
                Range = initZ + [fliplr(-Range) 0 Range];
            case 'log'
                Range = logspace(log10(6*AS.CoarseSS/nbsteps),log10(AS.CoarseR),nbsteps);
                Range = initZ + [fliplr(-Range) 0 Range];
            otherwise
                error(['Unknown range method: ' AS.method])
        end
        
        % Acquisition:
        [AS.Z_meas, AS.Sharpness_meas] = Zstack_sharpness(Range,AS.box);

        
        % The function we will have to minimize:
        funOPT = @(s) AS.reffunction(AS.Z_meas-s(1))' ...  S(1) = shift in z
            - (AS.Sharpness_meas.*s(2) ... % s(2) = scale
            + s(3)); % S(3) = overall shift in sharpness
        
        % I first get a rough estimate of where the minimum might be:
        % (extremely inefficient but at least it works)
        c= [];
        bs = -50:.1:50;
        if isempty(AS.scaletime), scale = 1; Sshift = 0; 
        else scale = AS.scaletime(end); Sshift = AS.Sharpshifttime(end);
        end
            
        for ind1 = bs
            c(end+1) = sum(funOPT([ind1,scale,Sshift]).^2);
        end
        [~,i] = min(c);
        
        
        % And then I run the optimization algorithm with this acquisition as the starting point for the optimization:
        [params, resnorm, ~, ~, ~] = lsqnonlin(funOPT,[bs(i),1,0],[bs(1),.3,-max(AS.Sharpness_truth)/10],[bs(end),1.5,max(AS.Sharpness_truth)/10]);
        
        AS.NewZ = AS.init_Z+params(1);
        
        
        
        % Update the time vectors:
        AS.exectimes(end+1) = toc(AS.starttic);
        AS.measZtime(end+1) = AS.NewZ;
        AS.shifttime(end+1) = params(1);
        AS.scaletime(end+1) = params(2);
        AS.Sharpshifttime(end+1) = params(3);
        AS.resnormtime(end+1) = resnorm;
        
        % And finally I run the smoothing algorithm:
        AS.decisiontime(end+1) = AS.init_Z + smoothlastNwithrejection(AS.shifttime,5,2);
        
        
        % Let's see those results:
        disp(['[Autofocus] Function minimization: Z Shift = ' num2str(params(1)) ', Scaling = ' num2str(params(2)) ', Sharpness shift = ' num2str(params(3)) ', Normalized residual = ' num2str(resnorm) ])
        
        
        % Plot curve comparison:
        fa = findobj('type','figure','name','autofocus');
        if fa
            set(0, 'CurrentFigure', fa);
        else
           figure('name','autofocus');
        end
        cla
        title('Autofocus')
        hold on;
        plot(AS.Z_meas,AS.Sharpness_meas.*params(2)+params(3),'bx');
        plot(AS.Z_truth+params(1),AS.reffunction(AS.Z_truth),'r');
        line([AS.NewZ AS.NewZ],[ylim],'Color','g','LineStyle','-');
        line([initZ initZ],[ylim],'Color','k','LineStyle','--');
        xlabel('Z position (um)')
        ylabel('Sharpness (a.u.)')
        xlim([-4 84]);
        
        
        
        % Plot the results over time:
        fa = findobj('type','figure','name','autofocusintime');
        if fa
            set(0, 'CurrentFigure', fa);
        else
           figure('name','autofocusintime');
        end
        cla
        hold on
        title('Focus over time')
        plot(AS.exectimes./3600,AS.measZtime,'r.');
        plot(AS.exectimes./3600,AS.decisiontime,'b');
        xlabel('time (h)')
        ylabel('position (um)')
        
        
        % Set the piezo in position:
        mmc.setPosition('PIZStage',AS.NewZ);

    end
    guihs.autofocus.settings = AS;
    pause(.2) % ?

end

% SHarpness estimation fucnton
function [sharpness]=estimate_sharpness(G)
    sharpness = mean2(abs(imgradient(G)));
end

% Initial setup function
function [AS] = init_autofocus()
global mmc mmGUI
width=mmc.getImageWidth();
height=mmc.getImageHeight();
AS.setupratio = 4;
persistent default_settings


if isempty(default_settings) % If it is not empty, we do nothing, as it is a persistent variable.
    default_settings = {'1', '20', '150', 'log'};
end
% Launch preview window
mmGUI.enableLiveMode(true);

% Store initial Zs
% User sets focus with wheel, modifies settings if necessary and clicks on button when done
x = inputdlg({'step size','range (X microns)','size of ROI','method (linear or log)'},'Autofocus settings', 1, default_settings);
default_settings = x;
CoarseSS = str2double(x{1}); % Step size
CoarseR = str2double(x{2}); % Range
refS = str2num(x{3});
method = x{4};
SetupSS = CoarseSS/AS.setupratio;

mmGUI.enableLiveMode(false);

AS.init_Z = mmc.getPosition('PIZStage');
mmc.snapImage();
img = mmc.getImage();
img=reshape(img,[width,height]);


% Give the possibility to move the square around.
% Display the image and the rectangle
px = width/2;
py = height/2;
flag = 1;
while flag
    f = figure();
    set(f,'name','Select Cross')
    imshow(imadjust(img'));
    if flag == 1
        title('Select Cross on image for autofocus and drift (double click on rectangle when finished)')
    elseif flag == 2
        title('I said DOUBLE-CLICK ON RECTANGLE when finished!')
    end
    h = imrect(gca, [px-refS/2, py-refS/2, refS, refS]);
    setResizable(h,0);
    drawnow;
    pos = wait(h);
    if ishandle(f)
        close(f);
        drawnow;
    end
    if ~isempty(pos)
        flag = 0;
    else
        flag = 2;
    end 
end

% Now we've got the definitive coordinates
% Once we have the coordinate of center cross, we can define an image
% that will be used as reference:
AS.box = pos;
px = pos(1)+refS/2; 
py = pos(2)+refS/2;
AS.ref = img((px-refS/2):(px + refS/2), (py-refS/2):(py+refS/2));
AS.init_Sharpness = estimate_sharpness(double(AS.ref));

% Calculate the range:
nbsteps = (round((CoarseR/SetupSS)*1.5));
switch method
    case 'linear'
        Range = cumsum(repmat(SetupSS,1,(round((CoarseR/SetupSS)*1.5))));
        Range = AS.init_Z + [fliplr(-Range) 0 Range];
    case 'log'
        Range = logspace(log10(6*CoarseSS/nbsteps),log10(CoarseR),nbsteps);
        Range = AS.init_Z + [fliplr(-Range) 0 Range];
    otherwise
            error(['Unknown range method: ' method])
end

[AS.Z_truth, AS.Sharpness_truth] = Zstack_sharpness(Range,AS.box);
%  if there are duplicate valuies, get rid of them
[~, alluniind] = unique(AS.Z_truth);
duplicate_ind = setdiff(1:size(AS.Z_truth), alluniind);
AS.Z_truth = AS.Z_truth(setdiff(1:numel(AS.Z_truth), duplicate_ind));
AS.Sharpness_truth = AS.Sharpness_truth(setdiff(1:numel(AS.Sharpness_truth), duplicate_ind));

% return to initial Z
mmc.setPosition('PIZStage',AS.init_Z );

% Get a spline function estimate of the points:
AS.reffunction = fit(AS.Z_truth,AS.Sharpness_truth','smoothingspline');


% Show the sharpness curve
fa = findobj('type','figure','name','autofocus');
if fa
    set(0, 'CurrentFigure', fa);
else
   figure('name','autofocus');
end
cla
title('Autofocus')
plot(AS.reffunction,AS.Z_truth, AS.Sharpness_truth);
hold on
line([AS.init_Z AS.init_Z],[ylim],'Color','g');


AS.CoarseSS = CoarseSS;
AS.CoarseR = CoarseR;
AS.method = method;
AS.measZtime = [];
AS.errortime = [];
AS.exectimes = [];
AS.starttic = tic();
AS.shifttime = [];
AS.scaletime = [];
AS.Sharpshifttime = [];
AS.resnormtime = [];
AS.decisiontime = [];

AS.initialized = true;
end

% Z Acquisition function:
function [Zpos, Sharpness] = Zstack_sharpness(Zobj,box)
global mmc guihs
width=box(3);
height=box(4);
Zpos = zeros(numel(Zobj),1);

msg = ['Autofocus: Acquiring sharpness profile...\n[' repmat('#',1,0) repmat('-',1,50) ']' num2str(0) '/' num2str(numel(Zobj)) '\n'];
fprintf(msg);
mmc.setProperty('Camera-1','ClearMode','Clear Never') % Increase acquisition speed by stopping the clearing of the sensor during the stack
oldROI = mmc.getROI();
mmc.setROI(box(1),box(2),box(3),box(4)); % Make sure that there is no conversion to perform between hte box and the ROI???
% Open the shutter:
mmc.setAutoShutter(false);
mmc.setProperty('Light shutter (Vincent-D1)','Command','Open')
% Get in position:
mmc.setPosition('PIZStage',Zobj(1) - 0.1 );
mmc.waitForDevice('PIZStage');
% Start z-stack:
for indZ = 1:numel(Zobj);
    mmc.setPosition('PIZStage',Zobj(indZ));
    mmc.waitForDevice('PIZStage');
    Zpos(indZ) = mmc.getPosition('PIZStage');
    mmc.snapImage();
    imgC{indZ} = mmc.getImage(); 
    imgC{indZ} =reshape(imgC{indZ},[width,height]);
    Sharpness(indZ) = estimate_sharpness(double(imgC{indZ}));
    imgC{indZ} = [];
    fprintf(repmat('\b',1,numel(sprintf(msg))));
    msg = ['Autofocus: Acquiring sharpness profile...\n[' repmat('#',1,round(50*indZ/numel(Zobj))) repmat('-',1,50-round(50*indZ/numel(Zobj))) ']' num2str(indZ) '/' num2str(numel(Zobj)) '\n'];
    fprintf( msg );
end
mmc.setProperty('Light shutter (Vincent-D1)','Command','Close');
mmc.setAutoShutter(true);
mmc.setROI(oldROI.getX(),oldROI.getY(),oldROI.getWidth(),oldROI.getHeight());
mmc.setProperty('Camera-1','ClearMode','Clear Pre-Exposure')
end


% Utilisties
function output = smoothlastNwithrejection(input, span, thresholdinsigmas)

    if numel(input) < span
        span = numel(input);
    end

    smoothedin = mysmooth(input,span);
    sigma = std(input'-smoothedin); % Weird I know
    input(abs(input'-smoothedin) > thresholdinsigmas*sigma) = input(abs(input'-smoothedin) > thresholdinsigmas*sigma).*NaN;
    lastN = input((end+1-span):end);
    lastN = lastN(~isnan(lastN));
    output = mean(lastN);

    if isempty(output) || isnan(output)
         output = smoothedin(end);

    end
end

function output = mysmooth(input, span)

    output = smooth(input,span);

    for ind1 = 1:floor(span/2)
        output(ind1) = mean(input(1:(floor(span/2)+ind1-1)));
        output(end+1-ind1) = mean(input((end+1-(floor(span/2)+ind1)):end));
    end
end