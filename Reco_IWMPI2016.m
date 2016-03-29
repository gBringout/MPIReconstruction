% This is the script going with the abstract presented at the conference
% IWMPI 2016 in Lübeck
% see https://www.researchgate.net/publication/299099474_MPI_system_matrix_reconstruction_making_assumptions_on_the_imaging_device_rather_than_on_the_tracer_spatial_distribution?ev=prf_pub

%% 1. Loading the required external functions
disp('1. Load the required external functions');
clear all
close all

addpath(genpath(fullfile('.','regu')))

%% 2. Loading the data
disp('2. Load and Pre-process the data'); tic

% define the name of the data
path_files = fullfile('.','data','1');
filename_SM = fullfile(path_files,'SystemResponse.h5');
filename_SM_Empty = fullfile(path_files,'SystemResponse_empty.h5');
filename_Meas = fullfile(path_files,'Measurement.h5');

% chek if hte folder exist, other create it
if exist(path_files)==0
    mkdir(path_files)
end
% if the data are not present, download them
if exist(filename_Meas,'file')==0
    vMatlab = version;
    url = 'https://github.com/KsenijaGraefe/SingleSidedData/blob/master/Data/Measurement.h5?raw=true';
    if str2num(vMatlab(1:3))>=8.4
        %if the version of matlab is recent enough , directly save the files
        websave(filename_Meas,url)
    else
        %otherwise open the browser
        web(url)
        disp(sprintf('Please save the file in path %s',path_files))
    end
    % to obtain infos on the file, use the command: infoSM = h5info(filename_Meas);
end
if exist(filename_SM,'file')==0
    vMatlab = version;
    url = 'https://github.com/KsenijaGraefe/SingleSidedData/blob/master/Data/SystemResponse.h5?raw=true';
    if str2num(vMatlab(1:3))>=8.4
        %if the version of matlab is recent enough , directly save the files
        websave(filename_SM,url)
    else
        %otherwise open the browser
        web(url)
        disp(sprintf('Please save the file in %s',path_files))
    end
    % to obtain infos on the file, use the command: infoSM = h5info(filename_SM);
end
if exist(filename_SM_Empty,'file')==0
    vMatlab = version;
    url = 'https://github.com/KsenijaGraefe/SingleSidedData/blob/master/Data/SystemResponse_empty.h5?raw=true';
    if str2num(vMatlab(1:3))>=8.4
        %if the version of matlab is recent enough , directly save the files
        websave(filename_SM_Empty,url)
    else
        %otherwise open the browser
        web(url)
        disp(sprintf('Please save the file in path %s',path_files))
    end
    % to obtain infos on the file, use the command: infoSM = h5info(filename_SM_Empty);
end

% Read the S and the empty measurement
SM_particle =  hdf5read(filename_SM,'systemResponseFrequencies');
SM_empty =  hdf5read(filename_SM_Empty,'systemResponseFrequencies');

% Convert the data as complex values
S = ((SM_particle(1:2:end,:,:,:) + 1i*SM_particle(2:2:end,:,:,:)));
Empty = ((SM_empty(1:2:end,:,:,:) + 1i*SM_empty(2:2:end,:,:,:)));

size(S)
size(Empty)
clear 'SM_particle' 'SM_empty' 

%% 3. Verification
disp('3. Display the SM')

%display the system matrix
figure
for i=1:100
    subplot(10,10,i)
    imagesc(abs(S(:,:,1,i+200)))
    axis square
    colormap('gray')
end

%% 4. Set the parameter used for acquisition
disp('4. Set the parameter used for acquisition')

f_start = 25000*1.7;
deltaf = 789.1414;
FC_start = round(f_start/deltaf);

nbrTotalFC = 22178;
nbrTotalFCPerChannel = nbrTotalFC/2;
freq = (0:nbrTotalFCPerChannel-1)*deltaf;
nbrPixel1 = 15;
nbrPixel2 = 15;

%% 5. Shape the SM and seperate both channels
disp('5. Shape the SM and seperate both channels')

SM1 = reshape(S,[nbrPixel1*nbrPixel2 nbrTotalFC]);
SM1 = SM1(:,1:end/2);

SM2 = reshape(S,[nbrPixel1*nbrPixel2 nbrTotalFC]);
SM2 = SM2(:,end/2+1:end);

clear 'S'

%% 6. Calculation of the orthogonality map or gramian matrix (see abstract for references)
disp('6. Calculation of the orthogonality map')
disp('This step takes around 100 minutes and requires around 4GB of RAM, Please confirm by pressing a key')
pause
% this step took 100 minutes on a Intel Xeon E5-4657L v2
S = [SM1 SM2];
nbrPixel = size(S,1);
nbrTotalFC = size(S,2);

orthogonalityMap = zeros(nbrTotalFC,nbrTotalFC);

% normalization
for i=1:nbrTotalFC
    nSMn(:,i) = S(:,i)/norm(S(:,i));
end
% orthogonality map (is symetrical)
for i=1:nbrTotalFC
    for j=i:nbrTotalFC
        orthogonalityMap(i,j) = abs(dot(nSMn(:,i),nSMn(:,j)));
    end
end
% symetric part
for i=1:nbrTotalFC
    for j=1:i
        orthogonalityMap(i,j) = orthogonalityMap(j,i);
    end
end

% display some properties
figure;plot(sum(orthogonalityMap));title('Sum Ortho')
figure;plot(mean(orthogonalityMap));title('Mean Ortho')
figure;plot(std(orthogonalityMap));title('Std Ortho')

%% 7. calculation of the energy in the SM
disp('7. Calculation of the enrgy of the SM')

w1 = zeros(1,nbrTotalFCPerChannel);
w2 = zeros(1,nbrTotalFCPerChannel);
for j= 1:nbrTotalFCPerChannel
  w1(j) = sum(abs(SM1(:,j)));
  w2(j) = sum(abs(SM2(:,j)));
end

figure;semilogy(freq,w1); hold all; semilogy(freq,w2)

%% 8. Load the phantom
disp('8. Load a phantom')

% if the data are not present, download them
Signal_hdf5 =  hdf5read(filename_Meas,'frequencies');
measurementNumber = 6; %Fig 15 on the right from the TMI paper (doi:10.1109/TMI.2015.2507187)

% put back the real and imaginary parts together
Signal = ((Signal_hdf5(:,measurementNumber,1:2:end) + 1i*Signal_hdf5(:,measurementNumber,2:2:end)));

% Seperated both channels
Signal1 = squeeze(Signal(1:end/2)).';
Signal2 = squeeze(Signal(end/2+1:end)).';

% Create the equivalent phantom
phantom = zeros(nbrPixel1,nbrPixel2);
phantom([11,13,15],6) = 1;
phantom(15,10) = 1;

figure; imagesc(phantom); axis square; colormap('gray')

clear 'Signal_hdf5' 'Signal'

%% 9. Calculate the measerements SNR
disp('9. Calculate the measerements SNR')

noise_mean1 = zeros(1,nbrTotalFCPerChannel);
noise_std1 = zeros(1,nbrTotalFCPerChannel);
noise_mean2 = zeros(1,nbrTotalFCPerChannel);
noise_std2 = zeros(1,nbrTotalFCPerChannel);

for i=1:nbrTotalFCPerChannel %for each FC
    noise_measure1 = Empty(:,:,:,i);
    noise_measure2 = Empty(:,:,:,nbrTotalFCPerChannel+i);
    % we measured it 225 times
    noise_mean1(i) = abs(mean(noise_measure1(:)));
    noise_std1(i) = abs(std(noise_measure1(:)));
    noise_mean2(i) = abs(mean(noise_measure2(:)));
    noise_std2(i) = abs(std(noise_measure2(:)));
end
clear noise_measure1 noise_measure2

figure;semilogy(noise_std1); hold all; semilogy(noise_std2)
figure;semilogy(noise_mean1); hold all; semilogy(noise_mean2)


Signal1_SNR_s_4 = (abs(Signal1)-noise_mean1)./noise_std1;
Signal1_SNR_s_4(1:FC_start) = 0;
Signal2_SNR_s_4 = (abs(Signal2)-noise_mean2)./noise_std2;
Signal2_SNR_s_4(1:FC_start) = 0;

figure;plot(freq/1000,Signal1_SNR_s_4); hold all; plot(freq/1000,Signal2_SNR_s_4)

%% 10. Reco with selected FC
disp('9. Reco with the presented method')

% Apply the thresdols on the orthogonality
idx_ortho = std(orthogonalityMap)<0.06 & mean(orthogonalityMap)<0.07; sum(idx_ortho)

% Apply the thresdols on the energy of the SM
idx1_energy = w1>0.010;sum(idx1_energy)
idx2_energy = w2>0.010;sum(idx2_energy)

% Apply the thresdols on the signal SNR
idx1_SNR = Signal1_SNR_s_4>10; sum(idx1_SNR)
idx2_SNR = Signal2_SNR_s_4>10; sum(idx2_SNR)

% total used frequency components
idx_reco = idx_ortho & [idx1_SNR idx2_SNR] & [idx1_energy idx2_energy]; sum(idx_reco)
%idx_reco = [idx1_SNR idx2_SNR] & [idx1_energy idx2_energy]; sum(idx_reco)

% form the truncated SM and signal
tSM = [SM1(:,idx_reco(1:nbrTotalFCPerChannel)) SM2(:,idx_reco(nbrTotalFCPerChannel+1:end))].';
tSignal = [Signal1(idx_reco(1:nbrTotalFCPerChannel)) Signal2(idx_reco(nbrTotalFCPerChannel+1:end))].';

% reconstruct
maxIterationReco = 50;
tic
C = artGael(tSM,tSignal,maxIterationReco);
toc

clear error1 error2
figure
subplot(1,2,1)
for i=1:maxIterationReco
    res = reshape(C(:,i),[nbrPixel1,nbrPixel2]);
    if i==1
        oldRes = res;
    end
    
    subplot(1,2,1)
    imagesc(real(res));
    
    error1(i) = sqrt((mean((real(oldRes(:))-real(res(:))).^2)));
    error2(i) = sqrt((mean((real(oldRes(:))-real(res(:))).^2)));
    axis square
    colormap('gray')
    %pause(1/25)
    oldRes = res;
end
subplot(1,2,2)
plot(error2);