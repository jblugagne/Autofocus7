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
        
        % TODO: Add XY correction to this function!!!!! (x-corr with the
        % reference.) (I could do an XY-drift correction function that simply calls this function.)
        
        % TODO: A new interface for homemade autofocus scripts.

        % If it has, acquire sharpness curve: 
        initZ = mmc.getPosition('PIZStage');    
        tic
        pause(.1)
        Range = cumsum(repmat(AS.CoarseSS,1,(round(.5*AS.CoarseR/AS.CoarseSS))));
        Range = initZ + [fliplr(-Range) 0 Range];
        [AS.Z_meas, AS.Sharpness_meas] = Zstack_sharpness(Range,AS.box);
        toc    

        % Cross correlate with the groundtruth: 
            % We have to zero-pad before (step size of groundtruth is ~3 times smaller than that of measurement plus we do the groundtruth over a range 3 times larger)
%         meas = zeros(1,AS.setupratio*numel(AS.Sharpness_meas));
%         meas(1:AS.setupratio:end) = AS.Sharpness_meas(:);
%         c = xcorr(AS.Sharpness_truth),meas);
%         c = c(numel(c)-(numel(AS.Sharpness_truth)+numel(meas)-2):end); % Matlab gives too many elements, I don't want most of them.
%         c = c(numel(meas):end-(numel(meas)-1));
        

        meas = AS.Sharpness_meas;
        for ind1 = 1:(numel(AS.Sharpness_truth)-AS.setupratio*numel(meas))
            c(ind1) = sum((meas-AS.Sharpness_truth(ind1:AS.setupratio:AS.setupratio*(numel(meas)-1)+ind1)).^2);
        end
        % Find max of Xcorrelation and compute position of user-defined
        % focus
 
        [~,imax] = min(c);
        truthoffset = AS.Z_meas(1) - AS.Z_truth(imax);
        AS.NewZ = AS.init_Z+truthoffset;

        fa = findobj('type','figure','name','autofocus');
        if fa
            set(0, 'CurrentFigure', fa);
        else
           figure('name','autofocus');
        end
        cla
        title('Autofocus')
        plot(AS.Z_truth,AS.Sharpness_truth,'m--');
        hold on;
        plot(AS.Z_meas,AS.Sharpness_meas,'bx');
        plot(AS.Z_truth+truthoffset,AS.Sharpness_truth,'m');
        line([AS.NewZ AS.NewZ],[ylim],'Color','g');
        line([initZ initZ],[ylim],'Color','r');
        
        
        AS.exectimes(end+1) = toc(AS.starttic);
        AS.measZtime(end+1) = AS.NewZ;
        AS.errortime(end+1) = AS.NewZ - initZ;
        if numel(AS.exectimes) > 1
            smootherror = smooth(AS.exectimes, AS.errortime, 1, 'rloess');
            AS.setZtime(end+1) = smootherror(end) + initZ;
        else
            AS.setZtime(end+1) = AS.NewZ;
        end
        
        fa = findobj('type','figure','name','autofocusintime');
        if fa
            set(0, 'CurrentFigure', fa);
        else
           figure('name','autofocusintime');
        end
        cla
        hold on
        title('Focus over time')
        plot(AS.exectimes./3600,AS.measZtime,'x');
        plot(AS.exectimes./3600,AS.measZtime-AS.errortime,'r.');
        plot(AS.exectimes./3600,AS.setZtime,'g--');
        
        xlabel('time (h)')
        ylabel('position (um)')
        
        
        mmc.setPosition('PIZStage',AS.NewZ);

    end
    guihs.autofocus.settings = AS;
    pause(.2)

end

% SHarpness estimation fucnton
function [sharpness]=estimate_sharpness(G)
    sharpness = mean2(abs(imgradient(G)));
end

% Initial setup function
function [AS] = init_autofocus()
global mmc mmGUI guihs
width=mmc.getImageWidth();
height=mmc.getImageHeight();
AS.setupratio = 4;
persistent default_settings


if isempty(default_settings) % If it is not empty, we do nothing, as it is a persistent variable.
    default_settings = {'1', '20', '150'};
end
% Launch preview window
mmGUI.enableLiveMode(true);

% Store initial Zs
% User sets focus with wheel, modifies settings if necessary and clicks on button when done
x = inputdlg({'step size','range (X microns)','size of ROI'},'Autofocus settings', 1, default_settings);
default_settings = x;
CoarseSS = str2double(x{1}); % Step size
CoarseR = str2double(x{2}); % Range
refS = str2num(x{3});
SetupSS = CoarseSS/AS.setupratio;

mmGUI.enableLiveMode(false);

AS.init_Z = mmc.getPosition('PIZStage');
% % < THE DIRTIEST HACK IN HISTORY >
% if strcmp('Light shutter (Vincent-D1)', mmc.getShutterDevice), fprintf(guihs.dirtylamphack,'OPEN'); end
% % </ THE DIRTIEST HACK IN HISTORY >
mmc.snapImage();
% % < THE DIRTIEST HACK IN HISTORY >
% fprintf(guihs.dirtylamphack,'CLOZ')
% % </ THE DIRTIEST HACK IN HISTORY >
img = mmc.getImage();
img=reshape(img,[width,height]);



% Crosses detection: I may redo that later. For now, it will be user
% selection
% R = imread('C:\Users\lab 513\JB\Autofocus5_Xfocus\CrossRef.tif');
% Re = edge(R,'log'); % Detect the edges (so we are not sensitive to lighting conditions)
% Ie = edge(img,'log');
% Rem = imclose(Re,strel('Disk',10)); % Do a math morpho closing to get a wider/more general form (Actually I could just replace it with a simple cross image, white on black...)
% Iem = imclose(Ie,strel('Disk',10));
% xc = xcorr2(double(Iem),double(Rem)); % xcorrelate the images to find the crosses
% There should be between 4 and 5 peaks, aligned and distributed equally:





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

    
% go +-range*1.5 en finestepsize
Range = cumsum(repmat(SetupSS,1,(round((CoarseR/SetupSS)*1.5))));
Range = AS.init_Z + [fliplr(-Range) 0 Range];
[AS.Z_truth, AS.Sharpness_truth] = Zstack_sharpness(Range,AS.box);

% return to initial Z
mmc.setPosition('PIZStage',AS.init_Z );

% Show the sharpness curve
fa = findobj('type','figure','name','autofocus');
if fa
    set(0, 'CurrentFigure', fa);
else
   figure('name','autofocus');
end
cla
title('Autofocus')
plot(AS.Z_truth, AS.Sharpness_truth,'m-');
hold on
plot(AS.init_Z,AS.init_Sharpness,'mx');


AS.CoarseSS = CoarseSS;
AS.CoarseR = CoarseR;
AS.initialized = true;
AS.measZtime = [];
AS.setZtime = [];
AS.errortime = [];
AS.exectimes = [];
AS.starttic = tic();
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
% % < THE DIRTIEST HACK IN HISTORY >
% fprintf(guihs.dirtylamphack,'OPEN');
% % </ THE DIRTIEST HACK IN HISTORY >
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
% % < THE DIRTIEST HACK IN HISTORY >
% fprintf(guihs.dirtylamphack,'CLOZ');
% % </ THE DIRTIEST HACK IN HISTORY >
mmc.setAutoShutter(true);
mmc.setROI(oldROI.getX(),oldROI.getY(),oldROI.getWidth(),oldROI.getHeight());
mmc.setProperty('Camera-1','ClearMode','Clear Pre-Exposure')
end