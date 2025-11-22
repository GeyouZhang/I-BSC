% This is the script file for motion-error-free 3D reconstruction according to the Binomial Self-compensation (BSC) technique. 
% Line 12: Switch between different sets of data
% Author: Geyou Zhang, University of Electronic Science and Technology of
% China, 2025/11/22
clc; clear all; close all;
%% Parameters Setting
addpath('./Package/');
load('../../Data/mCamera1Rectified.mat');
load('../../Data/mCamera2Rectified.mat');
load('../../Data/mProjector.mat');
% Select data
sFolderL = '../../Data/Statue/1_Rectified/';
sFolderR = '../../Data/Statue/2_Rectified/';
% sFolderL = '../../Data/Tissue/1_Rectified/';
% sFolderR = '../../Data/Tissue/2_Rectified/';
NSet = 4; FSet = 912/36;
iCameraWidth = 640; iCameraHeight = 480;
iFrameTotal = 50;
% Set Depth Range:
Zmin = -150; Zmax = 0;
% Set binomial order of motion error compensation, binomial order = 0 represents no compensation.
iBinomialOrder = 2;
iImageNum = iBinomialOrder + 4;
iBcThresh = 5;
% Compute Disparity Range
[ mDispMin, mDispMax ] = Func_DispartiyRange( Zmin, Zmax, mCamera1Rectified, mCamera2Rectified, iCameraHeight, iCameraWidth );

%% Image Sequence Loading
vmIL = nan( iCameraHeight, iCameraWidth, iFrameTotal );
vmIR = nan( iCameraHeight, iCameraWidth, iFrameTotal );
for i = 1:iFrameTotal  
    vmIL(:,:,i) = double( imread( sprintf( '%s%04d.bmp', sFolderL, i - 1 ) ) );
    vmIR(:,:,i) = double( imread( sprintf( '%s%04d.bmp', sFolderR, i - 1 ) ) );
end
%% Binomial Self-Compensation VS Traditional Four-step phase shifting for Dynamic 3D Scanning
% The computation speed can be accelerated by setting a larger number of workers in the MATLAB parallel pool
fig = figure;
set(gcf, 'Position', [0 0 1500 600]); %tiledlayout(2, 2, 'TileSpacing', 'none', 'Padding', 'none');
for i = 1:iFrameTotal - iImageNum + 1
    %% 3D reconstruction with traditional four-step phase shifting
    [ mPhaseWrapLeftFourStep ] = Func_PBSC( vmIL(:, :, i:i + 3), iBcThresh );
    [ mPhaseWrapRightFourStep ] = Func_PBSC( vmIR(:, :, i:i + 3), iBcThresh );
    % Correct Inherent Phase Shift
    dOffset = pi/2*mod(i - 1,4);
    mPhaseWrapLeftFourStep = mod(mPhaseWrapLeftFourStep + dOffset,2*pi);
    mPhaseWrapRightFourStep = mod(mPhaseWrapRightFourStep + dOffset,2*pi);    
    % Stereo Phase Unwarpping and 3D Reconstruction
    [ mXFourStep, mYFourStep, mZFourStep] = Func_Compute3D_PSPU( mPhaseWrapLeftFourStep, mPhaseWrapRightFourStep, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, 0 );
    % Remove the Outliers
    [ mZFourStep ] = Func_FiltZOutlier( mZFourStep, 9, 2 );   
    %% 3D reconstruction with P-BSC
    % Test PBSC
    [ mPhaseWrapLeft_PBSC ] = Func_PBSC( vmIL(:, :, i:i + iImageNum - 1), iBcThresh );
    [ mPhaseWrapRight_PBSC ] = Func_PBSC( vmIR(:, :, i:i + iImageNum - 1), iBcThresh );
    % Correct Inherent Phase Shift
    dOffset = pi/2*mod(i - 1,4);
    mPhaseWrapLeft_PBSC = mod(mPhaseWrapLeft_PBSC + dOffset,2*pi);
    mPhaseWrapRight_PBSC = mod(mPhaseWrapRight_PBSC + dOffset,2*pi);
    % Paraxial Stereo Phase Unwarpping and 3D Reconstruction
    [ mX_PBSC, mY_PBSC, mZ_PBSC ] = Func_Compute3D_PSPU( mPhaseWrapLeft_PBSC, mPhaseWrapRight_PBSC, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, 0 );
    % Remove the Outliers
    [ mZ_PBSC ] = Func_FiltZOutlier( mZ_PBSC, 9, 2 );   
    %% 3D reconstruction with I-BSC
    % Test IBSC
    [ mPhaseWrapLeft_IBSC ] = Func_IBSC( vmIL(:, :, i:i+ iImageNum - 1), iBcThresh );
    [ mPhaseWrapRight_IBSC ] = Func_IBSC( vmIR(:, :, i:i+ iImageNum - 1), iBcThresh );
    % Correct Inherent Phase Shift
    dOffset = pi/2*mod(i - 1,4);
    mPhaseWrapLeft_IBSC = mod(mPhaseWrapLeft_IBSC + dOffset,2*pi);
    mPhaseWrapRight_IBSC = mod(mPhaseWrapRight_IBSC + dOffset,2*pi);    
    % Paraxial Stereo Phase Unwarpping and 3D Reconstruction
    [ mX_IBSC, mY_IBSC, mZ_IBSC ] = Func_Compute3D_PSPU( mPhaseWrapLeft_IBSC, mPhaseWrapRight_IBSC, mDispMin, mDispMax, mCamera1Rectified, mCamera2Rectified, mProjector, FSet, 0 );
    % Remove the Outliers
    [ mZ_IBSC ] = Func_FiltZOutlier( mZ_IBSC, 9, 2 );      
    
    %% Draw the point clouds   
    disp(['Frame no.',num2str(i), '-----------------------------------------------------------------------------------------------------------']);       
    % Compute the depth range for visualization
    mu = mean( mZFourStep(:), 'omitnan' ); sigma = std( mZFourStep(:), 'omitnan' ); CMin = mu - 2*sigma; CMax = mu + 2*sigma;
      
    sgtitle(['Frame no.',num2str(i)], 'FontSize', 20, 'FontWeight', 'bold', 'Color',[1,1,1]);
    subplot(131);
    pcshow( [mXFourStep(:),mYFourStep(:),mZFourStep(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]); colormap(jet);title('Traditional four-step', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]); 
    
    subplot(132);
    pcshow( [mX_PBSC(:),mY_PBSC(:),mZ_PBSC(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]);colormap(jet); title('P-BSC', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]);    
    
    subplot(133);
    pcshow( [mX_IBSC(:),mY_IBSC(:),mZ_IBSC(:)] ); axis image; zlim([Zmin, Zmax]); xlim([-100, 40]); ylim([-120, 20]); 
    view([0 -90]);colormap(jet); title('I-BSC', 'FontSize', 20, 'FontWeight', 'bold'); caxis([CMin, CMax]);    
    
    h = axes(fig,'visible','off'); 
    cb = colorbar(h,'Position',[0.92 0.168 0.016 0.7], 'FontSize', 20, 'FontWeight', 'bold'); 
    set(get(cb,'Title'),'string','mm','FontWeight', 'bold', 'FontSize', 20, 'Color', 'white'); set(cb, 'FontWeight', 'bold', 'FontSize', 20, 'Color', 'white'); 
    caxis(h, [CMin, CMax]);
    drawnow;
end