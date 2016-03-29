%% 1. Loading the required external functions
disp('1. Load the required external functions');
clear all
close all

addpath(genpath(fullfile('.','regu')))
addpath(genpath(fullfile('.','ourFunctions')))

%% 2. Loading the data
disp('2. Load the data'); tic

% define the name of the data
filename_SM = fullfile('.','data','2','systemMatrix.h5');
filename_Meas = fullfile('.','data','2','measurement.h5');

% if the data are not present, download them
if exist(filename_SM,'file')==0
    urlwrite('http://media.tuhh.de/ibi/mdf/systemMatrix.h5',filename_SM)
    % to obtain infos on the file, use the command: infoSM = h5info(filename_SM);
    % or read the format documentation available at https://github.com/MagneticParticleImaging/MDF
end

if exist(filename_Meas,'file')==0
    urlwrite('http://media.tuhh.de/ibi/mdf/measurement_5.h5',filename_Meas)
    % to obtain infos on the file, use the command: infoSM = h5info(filename_Meas);
    % or read the format documentation available at https://github.com/MagneticParticleImaging/MDF
end

% For the System matrix (later named SM)
% read the data, saved as real numbers
S = h5read(filename_SM, '/calibration/dataFD');

% reinterpret as complex numbers
S = squeeze(S(1,:,:,:) + 1i*S(2,:,:,:));

% For the measurements
% the filename
% read and convert the data as complex numbers
% note that these data contain 500 measurements as of 20.03.2016
u = h5read(filename_Meas, '/measurement/dataFD');
u = squeeze(u(1,:,:,:) + 1i*u(2,:,:,:));
toc

%% 3. Pre-process and display the SM
disp('3. Pre-process and display the SM'); tic

% read the number of frequencies per channel for the SM
freq_SM = h5read(filename_SM, '/acquisition/receiver/frequencies');
numberFreq_SM = size(freq_SM,1);

% read the numbers of points used to discretize the 3D volume
number_Position = h5read(filename_SM, '/calibration/size');

% display one part of the first channel of the SM
figure
for i=1:100
    subplot(10,10,i)
    frequencyComponent = 50+i;
    imagesc(reshape(abs(S(:,frequencyComponent,1)),number_Position(1),number_Position(2)));
    axis square
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    title(sprintf('%i FC',frequencyComponent));
end
colormap(gray)
toc

%% 4. Pre-process and display a measurement
disp('4. Pre-process and display the SM'); tic

% load the frequencies corresponding to the matrices indexes
% hoping that they are the same for the SM and the phantom measurements :)
freq_Meas = h5read(filename_Meas, '/acquisition/receiver/frequencies');
numberFreq_Meas = size(freq_Meas,1);

figure
semilogy(freq_Meas,abs(u(:,1,1)))
title('Absolute value of a transformed FFT of the first measure on the first channel')
ylabel('Transformed FFT (unknown unit)')
xlabel('Frequency (Hz)')
toc

%% 5. Remove the frequencies which are lower than 30 kHz, as they are unreliable due to the anologue filter in the scanner
disp('5. Post-processing: remove the frequencies'); tic

% we supose that the same frequencies are measured on all channel for 
% the SM and the measurements
idxFreq = freq_Meas > 30e3;
S_truncated = S(:,idxFreq,:);
u_truncated = u(idxFreq,:,:);
toc

%% 6. Averaged the measurement used for the reconstruction over all temporal frames
disp('6. Post-processing: average the measurements'); tic

u_mean_truncated = mean(u_truncated,3);

%% 7. Make five simple reconstructions
disp('7. Make 5 simple recontruction');

% with the build in least square
% using a maximum of 1000 iterations
tic
maxIteration = 1000;
% and a small tolerance
tolerance = 10^-6;
c_lsqr = lsqr(S_truncated(:,:,1).', u_mean_truncated(:,1),tolerance,maxIteration);
disp('Least square')
toc

% and an external ART function
% using a maximum of 3 iterations
tic
maxIteration = 3;
c_art = art(S_truncated(:,:,1).',u_mean_truncated(:,1),maxIteration);
disp('ART')
toc

% and a modified version of the external ART function
% forcing a real and non-negative solution
% using a maximum of 3 iterations
tic
maxIteration = 3;
c_artGael = artGael(S_truncated(:,:,1).',u_mean_truncated(:,1),maxIteration);
disp('Modified ART')
toc

% and a normalized regularized kaczmarz approach
tic
maxIteration = 1;
c_normReguArt = regularizedKaczmarz(S_truncated(:,:),...
                        u_mean_truncated(:),...
                        maxIteration,...
                        1*10^-6,0,1,1);% lambda,shuffle,enforceReal,enforcePositive
disp('Regularized ART')
toc
                    
% and an regularized pseudoinverse approach
tic
[U,Sigma,V] = csvd(S_truncated(:,:).');
lambd = 5*10^3;
c_pseudoInverse = regularizedPseudoinverse(U,Sigma,V,u_mean_truncated,lambd,1,1);
disp('Pseudoinverse')
toc
%% 8. Display an image
disp('8. Display the 5 reconstruction')

figure
subplot(3,2,1)
imagesc(real(reshape(c_lsqr(:),number_Position(1),number_Position(2))));
colormap(gray); axis square
title({'Matlab least square - 1 channel';'1000th iterations / real part'})

subplot(3,2,2)
imagesc(real(reshape(c_art(:,1),number_Position(1),number_Position(2))));
colormap(gray); axis square
title({'External ART - 1 channel';'1st iterations / real part'})

subplot(3,2,3)
imagesc(real(reshape(c_artGael(:,1),number_Position(1),number_Position(2))));
colormap(gray); axis square
title({'Modified ART - 1 channel';'1st iterations / real part'})

subplot(3,2,4)
imagesc(real(reshape(c_normReguArt(:),number_Position(1),number_Position(2))));
colormap(gray); axis square
title({'Regularized and modified ART - 3 channels';'1 iterations / lambda = 10^{-6} / real part'})

subplot(3,2,5)
imagesc(real(reshape(c_pseudoInverse(:),number_Position(1),number_Position(2))));
colormap(gray); axis square
title({'Pseudoinverse - 3 channels';' lambda = 5*10^{3} / real part'})
