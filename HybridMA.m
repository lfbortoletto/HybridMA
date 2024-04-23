% This function corrects motion artifacts in time series data using spline 
% interpolation. The function identifies segments where motion artifacts 
% exceed a specified threshold and applies a spline-based correction.
%
% The algorithm follows the procedure found in https://doi.org/10.1117/1.NPh.7.1.015001.
%
% Inputs:
%   timeSeriesData: Matrix of input data (time points x channels)
%   acquisitionFrequency: Sampling frequency of the data acquisition
%       system.
%   threshold: Motion detection threshold.
%   K: Free parameter for motion correction (leave empty for default. Default: 2.5*acquisitionFrequency).
%
% Outputs:
%   hybridCorrectedData: Motion-corrected time series by hybrid approach, using spline and wavelet 
%       combination.
%
% Created by: S. L. Novi (2017/08/01)
%
% Last Modified by:
%   L. F. Bortoletto (2024/04/18): code structure optimization.
%   L. F. Bortoletto (2024/04/23): add wavelet to motion correction.
% -------------------------------------------------------------------------

function hybridCorrectedData = HybridMA(timeSeriesData, acquisitionFrequency, splineThreshold, waveletThreshold, K)
    if isempty(K)
        % Use 'K' default value.
        K = round(2.5 * acquisitionFrequency);
    end
    % Correct data using spline interpolation.
    splineCorrectedData = MotionCorrectionMARA(timeSeriesData, acquisitionFrequency, splineThreshold, K);
    
    % Correct data using wavelet transform.
    SD.MeasListAct = ones(size(timeSeriesData, 2), 1); % List of channels to be corrected.
    hybridCorrectedData = hmrMotionCorrectWavelet(splineCorrectedData, SD, waveletThreshold);
end

function correctedTimeSeries = MotionCorrectionMARA(timeSeriesData, acquisitionFrequency, motionThreshold, K)

    windowSize = (2 * K) + 1;
    numTimePoints = size(timeSeriesData, 1);
    movingStdDev = zeros(numTimePoints, size(timeSeriesData, 2));
    
    for chan = 1:size(timeSeriesData, 2)
        for t = K + 1:numTimePoints - K
            localSegment = timeSeriesData(t-K : t+K, chan);
            % Calculate moving standard deviation
            termA = localSegment' * localSegment;
            termB = (sum(localSegment))^2 / windowSize;
            movingStdDev(t, chan) = sqrt((termA - termB) / windowSize);
        end
        
        % Threshold based on standard deviation criteria
        meanStd = mean(movingStdDev(K+1:end, chan));
        stdStd = std(movingStdDev(K+1:end, chan));
        movingStdDev(movingStdDev(:, chan) < meanStd + motionThreshold * stdStd, chan) = 0;
    end
    
    % Identify segments with motion artifacts
    motionArtifactSegments = identifyMotionSegments(movingStdDev);
    
    % Reconstruct time series correcting the artifacts
    correctedTimeSeries = reconstructTimeSeries(timeSeriesData, motionArtifactSegments, acquisitionFrequency);
end

function motionArtifactSegments = identifyMotionSegments(movingStdDev)
    % Identify continuous segments of the time series with motion artifacts.
    % Inputs:
    %   movingStdDev: Matrix of moving standard deviations for each sample and channel.
    % Outputs:
    %   motionArtifactSegments: Cell array containing start and end indexes of motion artifacts per channel.
    
    motionArtifactSegments = cell(1, size(movingStdDev, 2));
    for chan = 1:size(movingStdDev, 2)
        activePoints = find(movingStdDev(:, chan) > 0);
        flippedActivePoints = flip(activePoints);
        if ~isempty(activePoints)
            segmentsBeggining = find(diff([0; activePoints; max(activePoints)+1]) > 1); % Motion artifact beggining indexes.
            segmentsBeggining = reshape(segmentsBeggining, numel(segmentsBeggining),[]);
            segmentsEnding = find(diff(flip([0; activePoints; max(activePoints)+1])) < -1); % Motion artifact ending indexes.
            segmentsEnding = reshape(segmentsEnding(1:end-1), numel(segmentsEnding(1:end-1)),[]);
            segments = sort([activePoints(segmentsBeggining); flippedActivePoints(segmentsEnding)]);
            motionArtifactSegments{chan} = reshape(segments, numel(segments),[]);
            motionArtifactSegments{chan}(end+1) = activePoints(end);
        end
    end
end

function correctedTimeSeries = reconstructTimeSeries(originalTimeSeries, segments, acquisitionFrequency)
    % Reconstructs the original time series by correcting identified motion artifacts.
    % Inputs:
    %   originalTimeSeries: Original matrix of time series data (time points x channels).
    %   segments: Cell array of segments indicating motion artifacts.
    %   acquisitionFrequency: Sampling frequency of the data acquisition
    %   system.
    % Output:
    %   correctedTimeSeries: Corrected time series data.
    
    MotionArtifactAdjusted = originalTimeSeries;
    timings = linspace(0, size(originalTimeSeries,1)*acquisitionFrequency, size(originalTimeSeries,1)); %% Create vector t.
    timings = timings';
    for chan = 1:size(originalTimeSeries, 2)
        if isempty(segments{chan})
            continue;
        end
        for i = 1:size(segments{chan}, 1) - 1
            % Spline interpolation over the segments of motion artifacts
            segmentIndices = segments{chan}(i):segments{chan}(i+1);
            if length(segmentIndices) <= 3
                % Zero out if segment is too small.
                MotionArtifactAdjusted(segmentIndices, chan) = 0;
            else
                % Perform spline interpolation and adjust the segment.
                pp = 0.01;
                splineTimeSeries = csaps(timings(segmentIndices),originalTimeSeries(segmentIndices, chan), pp, timings(segmentIndices));
                MotionArtifactAdjusted(segmentIndices, chan) = originalTimeSeries(segmentIndices, chan) - splineTimeSeries; 
            end
        end
    end
    
    % Correcting for the different baselines
    correctedTimeSeries = originalTimeSeries;
    for chan = 1:size(originalTimeSeries,2)
        for nSeg = 1:numel(segments{chan})
            if nSeg==1
                t1 = 1:segments{chan}(nSeg); %% First Segment (Frames).
                t2 = segments{chan}(nSeg):segments{chan}(nSeg+1); % Second Segment (Frames).
                x1 = originalTimeSeries(t1, chan);
                x2 = MotionArtifactAdjusted(t2, chan);               
            elseif mod(nSeg,2)==0 && nSeg~= numel(segments{chan})
                t1 = segments{chan}(nSeg-1):segments{chan}(nSeg); % p segment.
                t2 = segments{chan}(nSeg):segments{chan}(nSeg+1); % p+1 segment.
                x1 = correctedTimeSeries(t1,chan);
                x2 = originalTimeSeries(t2,chan);             
            elseif mod(nSeg,2) == 1
                t1 = segments{chan}(nSeg-1):segments{chan}(nSeg); % p segment.
                t2 = segments{chan}(nSeg):segments{chan}(nSeg+1); % p+1 segment.
                x1 = correctedTimeSeries(t1,chan);
                x2 = MotionArtifactAdjusted(t2,chan);
            elseif nSeg == numel(segments{chan})
                t1 = segments{chan}(nSeg-1):segments{chan}(nSeg); % p segment.
                t2 = segments{chan}(nSeg):size(originalTimeSeries,1);
                x1 = correctedTimeSeries(t1,chan);
                x2 = originalTimeSeries(t2,chan);
            end
            correctedTimeSeries(t2,chan) = CorrectBaselineMARA(x1, x2, acquisitionFrequency);
        end
    end
end

function [baselineAdjusted] = CorrectBaselineMARA(x1, x2, acquisitionFrequency)

%%% Criterias to correct the baseline of each time series segment. 
%%% All criterias were extract from doi:10.1088/0967-3334/31/5/004 
%%% in Table 1. 

a = round(acquisitionFrequency/3);
b = round(2*acquisitionFrequency);


if length(x1)<=a && length(x2)<=a
    v = mean(x1) - mean(x2);
   
elseif a<length(x1) && b>length(x1) && length(x2)<=a 
    
    v = mean(x1(end-a:end)) - mean(x2);
    %x3 = 
elseif length(x1)>=b && length(x2)<=a 
    l = length(x1);
    v = mean(x1(end-round(l/10):end)) - mean(x2);
    
elseif length(x1)<=a && length(x2)<b && length(x2)>a 
    v = mean(x1) - mean(x2(1:a));
    
elseif a<length(x1) && b>length(x1) && a<length(x2) && b>length(x2)
    v  = mean(x1(end-a:end)) - mean(x2(1:a));

elseif length(x1)>=b && length(x2)<b && length(x2)>a 
    l = length(x1);
    v = mean(x1(end-round(l/10):end)) - mean(x2(1:a));
    
elseif length(x1)<=a && length(x2)>=b
    l = length(x2);    
    v = mean(x1) - mean(x2(1:round(l/10)));
    
elseif a<length(x1) && b>length(x1) && length(x2)>=b
    l = length(x2);
    v = mean(x1(end-a:end)) - mean(x2(1:round(l/10)));
else
    l1 = length(x1);
    l2 = length(x2);
    v = mean(x1(end-round(l1/10):end)) - mean(x2(1:round(l2/10)));
end
    if isnan(v)
       v=0; 
    end
    baselineAdjusted = x2+v;
    
end

