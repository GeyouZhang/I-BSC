function [ mZ ] = Func_FiltZOutlier( mZ, iWindowSize, iThresh )
mNan = isnan(mZ);
mZ(mNan) = 0;
mZFilted = medfilt2(mZ, [iWindowSize,iWindowSize]);
vOutlier = find( abs(mZFilted - mZ) > iThresh);
mZ(vOutlier) = nan;
mZ(mNan) = nan;