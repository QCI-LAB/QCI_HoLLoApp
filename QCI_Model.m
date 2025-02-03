classdef QCI_Model < handle
    %QCI_MODEL is a bundle used for hologram creation and manipulation
    %   Detailed explanation goes here
    
    properties (Access = protected)
        %Initial Data
        Names                   (1,:) string = []
        Holograms               (:,:,:) double                   = []       % Hologram image (matrix)
        Wavelengths             (1,:) double {mustBeNonnegative} = []       % Wavelength of light used [micro meters]
        ZCoordinates            (1,:) double                     = []       % Distance from CCD to sample [micro meters] on the optical axis
        CameraPixelSize         (1,1) double = 1                            % Size of a pixel on the CCD camera used to obtain the hologram
        MediumRefractiveIndex   (1,1) double = 1                            % Refractive index of the medium in which the hologram was taken (assume air)

        % Derrived data
        YCoordinates            (1,:) double {mustBeNonnegative} = []       % Y hologram coordinates in physical units [micro meters]. Y axis points down from the left upper corner.
        XCoordinates            (1,:) double {mustBeNonnegative} = []       % X hologram coordinates in physical units [micro meters]. X axis points right from the left upper corner.
        YInverseCoordinates     (1,:) double                     = []       % Y hologram coordinates in inverse space;
        XInverseCoordinates     (1,:) double                     = []       % X hologram coordinates in inverse space;
        WaveNumbers             (1,:) double {mustBeNonnegative} = []       % Wavenumber in [1/micrometer = 1e6 m]
        FreeSpaceImpulseMatrixes(:,:,:) double                   = []
        HologramsFFT            (:,:,:)                          = []
    end
    
    methods
        function obj = QCI_Model(names, holograms, wavelengths, zCoordinates, pixelSize)
            if nargin == 0
                return
            end
            %QCI_MODEL Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Names = strings(1,length(names));
            obj.Names = names;
            
            for i = 1:length(wavelengths)
                obj.Holograms(:,:,i) = holograms(:,:,i);
            end

            obj.Wavelengths = wavelengths;
            obj.ZCoordinates = zCoordinates;
            obj.CameraPixelSize = pixelSize;
            [ySize, xSize] = size(obj.Holograms(:,:,1));
            
            % Deriving data
            obj.YCoordinates = (0:(ySize-1))*obj.CameraPixelSize;
            obj.XCoordinates = (0:(xSize-1))*obj.CameraPixelSize;
            inverseStepY = 1/ySize/pixelSize;
            inverseStepX = 1/xSize/pixelSize;
            obj.YInverseCoordinates =  -ySize/2*inverseStepY:inverseStepY:(ySize/2 - 1)*inverseStepY;
            obj.XInverseCoordinates =  -xSize/2*inverseStepX:inverseStepX:(xSize/2 - 1)*inverseStepX;
            obj.WaveNumbers = 2*pi./obj.Wavelengths;
            
            for i = 1:length(wavelengths)
                obj.HologramsFFT(:,:,i) = fft2(obj.Holograms(:,:,i));
                obj.FreeSpaceImpulseMatrixes(:,:,i) = fftshift(obj.WaveNumbers(i) *...
                    sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                    (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
            end
        end
        
        function names = getNames(obj, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(indexes == 0)
                names = obj.Names;
            else
                names = obj.Names(indexes);
            end
     
            return
        end

        function holograms = getHolograms(obj, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(indexes == 0)
                holograms = obj.Holograms;
            else
                holograms = obj.Holograms(:,:, indexes);
            end
        end

        function holograms = setHolograms(obj, holograms, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            if(size(holograms, 3) == size(obj.Holograms,3))
                    obj.Holograms = double(holograms);
                    obj.HologramsFFT = fft2(holograms);
            elseif(size(holograms, 3) == size(indexes))
                obj.Holograms(:,:,indexes) = double(holograms);
                obj.HologramsFFT(:,:,indexes) = fft2(holograms);
            else
                error("Dimension of indexes array does not match the amount of given holograms.");
            end
        end


        function addHologram(obj, name, hologram, wavelength, zCoordinate)
            obj.Names(end + 1) = name;
            obj.Holograms(:, :, end + 1) = hologram;
            obj.Wavelengths(end + 1) = wavelength;
            obj.ZCoordinates(end + 1) = zCoordinate;
            obj.WaveNumbers(end + 1) = 2*pi./obj.Wavelengths(end);
            
            obj.HologramsFFT(:,:,end + 1) = fft2(obj.Holograms(:,:,end));
            [ySize, xSize] = size(obj.Holograms(:, :, 1));
            obj.FreeSpaceImpulseMatrixes(:,:,end + 1) = fftshift(obj.WaveNumbers(end) *...
                sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(end)^2 *...
                (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
        end

        function removeHologram(obj, hologramIndex)
            obj.Holograms(:, :, hologramIndex) = [];
            obj.Wavelengths(hologramIndex) = [];
            obj.ZCoordinates(hologramIndex) = [];
            obj.WaveNumbers(hologramIndex) = [];
            
            obj.HologramsFFT(:, :, hologramIndex) = [];
            obj.FreeSpaceImpulseMatrixes(:, :, hologramIndex) = [];

            if(isempty(obj.Holograms))
                obj.XInverseCoordinates = [];
                obj.YInverseCoordinates = [];
                obj.XCoordinates = [];
                obj.YCoordinates = [];
                obj.CameraPixelSize = 1;
                obj.MediumRefractiveIndex = 1;
            end
        end

        function wavelengths = getWavelengths(obj, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(indexes == 0)
                wavelengths = obj.Wavelengths;
            else
                wavelengths = obj.Wavelengths(indexes);
            end
        end

        function setWavelengths(obj, wavelengths, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(size(wavelengths) == size(obj.Wavelengths))
                for i = 1:length(wavelengths)
                    if(wavelengths(i) ~= obj.Wavelengths(i))
                        obj.FreeSpaceImpulseMatrixes(:,:,i) = fftshift(obj.WaveNumbers(i) *...
                            sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                            (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
                    end
                end
                obj.Wavelengths = wavelengths;
            elseif(size(wavelengths) == size(indexes))  
                for i = indexes
                    if(wavelengths(i) ~= obj.Wavelengths(i))
                        obj.FreeSpaceImpulseMatrixes(:,:,i) = fftshift(obj.WaveNumbers(i) *...
                            sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                            (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
                    end
                end

                obj.Wavelengths(indexes) = wavelengths;
            else
                error("Dimension of indexes array does not match the amount of given wavelengths.");
            end

            return
        end

        function zCoordinates = getZCoordinates(obj, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(indexes == 0)
                zCoordinates = obj.ZCoordinates;;
            else
                zCoordinates = obj.ZCoordinates(indexes);
            end
        end

        function setZCoordinates(obj, zCoordinates, indexes)
            % Setter for ZCoordinates
            if(size(zCoordinates) == size(obj.ZCoordinates))
                obj.ZCoordinates = zCoordinates;
            elseif(size(zCoordinates) == size(indexes))  
                obj.ZCoordinates(indexes) = zCoordinates;
            else
                error("Dimension of indexes array does not match the amount of given zCoordinates.");
            end
        end
        
        % Operations on holograms
        function hologramPostPropagation = propagate(obj, hologramIndex, willReplace)
            % Propagates the image with Angular Scaling method
            kernelExponent = obj.FreeSpaceImpulseMatrixes(:,:,hologramIndex)*obj.ZCoordinates(hologramIndex);
            kernelExponent = kernelExponent - kernelExponent(1,1);
            %hologramPostPropagation = zeros(size(obj.Holograms,1), size(obj.Holograms,2));

            if(obj.ZCoordinates(hologramIndex)<0)
                kernel = exp(-1i*kernelExponent);
                
                FT_Uout = conj(obj.HologramsFFT(:,:, hologramIndex));
                FT_Uout = flip(flip(FT_Uout, 1), 2);  % Odbicie w pionie i poziomie
                FT_Uout = circshift(FT_Uout, [1, 1]); % Przesunięcie o 1 piksel

                propagatedFFT = kernel.*FT_Uout;

                hologramPostPropagation = conj(ifft2(propagatedFFT));
            else
                kernel = exp(1i*kernelExponent);
                propagatedFFT = kernel.*obj.HologramsFFT(:,:, hologramIndex);
                hologramPostPropagation = ifft2(propagatedFFT);
            end
            
            hologramPostPropagation = double(hologramPostPropagation);

            if(nargin == 3)
                if(willReplace)
                    obj.Holograms(:,:,hologramIndex) = hologramPostPropagation;
                    obj.HologramsFFT(:,:, hologramIndex) = propagatedFFT;
                end
            end
        end

        function [fixedReferenceImages, fixedHolograms, correctedHeights, tforms] = ThreePointShiftAndScalingCorrection(obj, referenceImages, points)
            
            for pointIndex = 1:3
                fixedPoints(pointIndex,:) = points{end}{pointIndex};
            end
            
            for tt = 1:size(referenceImages,3)-1
                for pointIndex = 1:3
                    movingPoints(pointIndex,:) = points{tt}{pointIndex};
                end

                tform = fitgeotrans(movingPoints,fixedPoints,'similarity');
                Roriginal = imref2d(size(referenceImages(:,:,end)));
                fixedHolograms(:,:,tt) = imwarp(obj.Holograms(:,:,tt),tform,'OutputView',Roriginal);
                fixedReferenceImages(:,:,tt) = imwarp(referenceImages(:,:,tt),tform,'OutputView',Roriginal);
                correctedHeights(tt) = round(obj.ZCoordinates(tt).*tform.T(1,1).^2);
                tforms{tt} = tform;
            end
            fixedHolograms(:,:,end+1) = obj.Holograms(:,:,end);
            fixedReferenceImages(:,:,end+1) = referenceImages(:,:,end);
            correctedHeights(end+1) = obj.ZCoordinates(end);
        end
        
        function [fixedReferenceImages, fixedHolograms, correctedHeights, tforms] = automaticalShiftAndScalingCorrection(obj, referenceImages)
        %   Align and scale holograms to reference images based on 
        %   SURF feature matching and similarity transformation.
        %
        %   Inputs:
        %       referenceImages  - 3D array of reference images (amplitude and phase).
        %       holograms        - 3D array of holograms to align with the references.
        %
        %   Outputs:
        %       fixedReferenceImages - Aligned and scaled reference images.
        %       fixedHolograms       - Aligned and scaled holograms.
        %       tforms               - Cell array of geometric transformations.
            % --- Initialization ---
            if(~isempty(obj.Holograms))
                holograms = obj.Holograms;
            else
                error("Module has no holograms to work with.");
            end

            imageNumber = size(referenceImages, 3);
            fixedReferenceImages = zeros(size(referenceImages), 'like', referenceImages);
            fixedHolograms = zeros(size(holograms), 'like', holograms);
            tforms = cell(1, imageNumber-1);  % Store transformations
            correctedHeights = zeros(1,imageNumber);
        
            % Select the last reference image as the alignment target
            referenceImageFixed = -angle(referenceImages(:, :, end));
            fixedReferenceImages(:, :, end) = referenceImageFixed;
            fixedHolograms(:, :, end) = holograms(:, :, end);
            correctedHeights(end) = obj.ZCoordinates(end);
            
            % Detect features in the fixed reference image
            featureDetectionThreshold = 10000;
            ptsFixed = detectSURFFeatures(referenceImageFixed, 'MetricThreshold', featureDetectionThreshold);
            [featuresFixed, validPtsFixed] = extractFeatures(referenceImageFixed, ptsFixed);
        
            % --- Loop through all images except the last one ---
            for k = 1:imageNumber-1
                % Current distorted reference image and hologram
                distortedImage = -angle(referenceImages(:, :, k));
                distortedHologram = holograms(:, :, k);
        
                % Step 1: Detect and extract features
                ptsDistorted = detectSURFFeatures(distortedImage, 'MetricThreshold', featureDetectionThreshold);
                [featuresDistorted, validPtsDistorted] = extractFeatures(distortedImage, ptsDistorted);
        
                % Step 2: Match features
                indexPairs = matchFeaturesInRadius(featuresFixed, featuresDistorted, ...
                                                   validPtsDistorted.Location, validPtsFixed.Location, 150);
        
                % Retrieve matched points
                matchedFixed = validPtsFixed(indexPairs(:,1));
                matchedDistorted = validPtsDistorted(indexPairs(:,2));
        
                % Step 3: Estimate geometric transformation
                [tform, inlierIdx] = estimateGeometricTransform2D(matchedDistorted, ...
                                                                 matchedFixed, 'similarity');
        
                % Step 4: Apply transformation to images and holograms
                Routput = imref2d(size(referenceImageFixed));
                fixedReferenceImages(:, :, k) = imwarp(distortedImage, tform, 'OutputView', Routput);
                fixedHolograms(:, :, k) = imwarp(distortedHologram, tform, 'OutputView', Routput);
                tforms{k} = tform;
                correctedHeights(k) = round(obj.ZCoordinates(k).*tform.T(1,1).^2);
        
                % Optional Debugging (uncomment for visualization)
                % figure; showMatchedFeatures(referenceImageFixed, distortedImage, ...
                %                            matchedFixed(inlierIdx), matchedDistorted(inlierIdx));
                % title(sprintf('Image %d Matching Points (Inliers Only)', k));
            end
        end

        function reconstruction = IGA(obj, iter, sigma)
        % Iterative Gabor Averagin (IGA) - method for phase retrieval from multiple
        % in-line holograms collected with different defocus and/or wavelength.
        % Algorithm was designed to work with low signal-to-noise-ratio data
        %
        % Inputs:
        %   H - 3D array containing registered in-line holograms (at least 2
        %       hologram are reguired)
        %   z - vector containing defocus distances for each hologram 
        %       (AS_propagate_p(H(:,:,n),z(n),lambda(n),dx) should give the
        %       in-focus reconstruction)
        %   lambda - vector containing wavelengths used to collect each in-line
        %       hologram
        %   dx - effective pixel size of the camera (camera pixel size divided by 
        %       system magnification)
        %   sigma - denoising factor that balances the twin image and shot noise
        %       removal. For larger sigma, shot noise should be minimized more
        %       effectively, while twin image is minimized better for smaller
        %       sigma. Default - sigma = 2; 
        %       sigma = 0 - GS method. sigma = inf - GA method. 
        %       Recommended values:
        %       sigma = 0 - noise-free data (only simulations)
        %       sigma = 1 - good quality data with insignificant shot noise
        %       sigma = 2 - regular or noisy data
        %       sigma = 4 - very strong shot noise, but twin image still present
        %       sigma = inf - shot noise larger than signal
        %   iter - number of iterations. Default - iter = 5;
        % Output:
        %   R - reconstructed complex optical field at the sample plane. 
        %       abs(R) -> amplitude; angle(R) -> phase.
        % 
        % More details at/cite as:
        %   M. Rogalski, P. Arcab, E. Wdowiak, J. Á. Picazo-Bueno, V. Micó, 
        %   M. Józwik, M. Trusiak, "Hybrid iterating-averaging low photon budget 
        %   Gabor holographic microscopy", Submitted 2024
        % 
        % Created by:
        %   Mikołaj Rogalski
        %   Warsaw University of Technology, Institute of Micromechanics and
        %   Photonics
        %   mikolaj.rogalski.dokt@pw.edu.pl
        % 
        % Last modified:
        %   04.09.2024
            
            if sigma == 0 % GS method
                reconstruction = obj.GerchbergSaxton(iter, sigma);
            elseif sigma == inf % GA method
                reconstruction = obj.GaborAveraging();
            else % IGA method
                % Combine GS and GA methods
                GerchbergSaxtonReconstruction = obj.GerchbergSaxton(iter, sigma);
                GaborAveragingReconstruction = obj.GaborAveraging();

                realPart = real(GerchbergSaxtonReconstruction) + real(GaborAveragingReconstruction) - imgaussfilt(real(GaborAveragingReconstruction), sigma);
                imaginaryPart = imag(GerchbergSaxtonReconstruction) + imag(GaborAveragingReconstruction) - imgaussfilt(imag(GaborAveragingReconstruction), sigma);
                reconstruction = realPart + 1i.*imaginaryPart;
            end
        end
        
        function reconstruction = GerchbergSaxton(obj, iter, sigma) 
            % Multi-height Gerchberg-Saxton in-line holography phase retrieval
            if(sigma > 0)
                hologramsUnblurred = obj.Holograms;
                hologramsBlurred = imgaussfilt(hologramsUnblurred, sigma);
                obj.Holograms = hologramsBlurred;
            end

            amplitudeImages = sqrt(obj.Holograms); % Amplitude = sqrt(intensity)
            hologramCount = size(obj.Holograms,3);
            isSingleWavelength = isscalar(unique(obj.Wavelengths));

            obj.addHologram("intermediateHologram", amplitudeImages(:, :, 1), obj.Wavelengths(1), obj.ZCoordinates(1)); % ZCoordinate will be replaced
            intermediateHologramIndex = hologramCount + 1;

            if(isSingleWavelength)
                for tt = 1:iter
                    for nn = 1:(hologramCount-1)

                        % Propagate to nn+1 hologram plane
                        propagationDistance = obj.ZCoordinates(nn+1) - obj.ZCoordinates(nn);
                        obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                        intermediateField = obj.propagate(intermediateHologramIndex);

                        % Actualize optical field with nn+1 amplitude
                        intermediateField = intermediateField./abs(intermediateField) .* amplitudeImages(:,:,nn+1);
                        obj.setHolograms(intermediateField, intermediateHologramIndex);
                    end
                    % Propagate to 1st hologram plane
                    propagationDistance = obj.ZCoordinates(1) - obj.ZCoordinates(hologramCount);
                    obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                    intermediateField = obj.propagate(intermediateHologramIndex);
                    % Actualize optical field with 1st amplitude
                    intermediateField = intermediateField./abs(intermediateField) .* amplitudeImages(:,:,1);
                    obj.setHolograms(intermediateField, intermediateHologramIndex);
                end
            else
                for tt = 1:iter
                    for nn = 1:(hologramCount-1)
                        % Propagate to object plane
                        propagationDistance = -obj.ZCoordinates(nn);
                        obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                        intermediateReconstruction = obj.propagate(intermediateHologramIndex);

                        % Rescale phase to nn+1 wavelength
                        phase = angle(intermediateReconstruction) * obj.Wavelengths(nn)/obj.Wavelengths(nn+1);
                        intermediateReconstruction = abs(intermediateReconstruction) .* exp(1i.*phase);
                        obj.setHolograms(intermediateReconstruction, intermediateHologramIndex);

                        % Propagate to nn+1 hologram plane
                        propagationDistance = obj.ZCoordinates(nn+1);
                        obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                        obj.setWavelengths(obj.Wavelengths(nn+1), intermediateHologramIndex);
                        intermediateField = obj.propagate(intermediateHologramIndex);

                        % Actualize optical field with nn+1 amplitude
                        intermediateField = intermediateField ./ abs(intermediateField) .* amplitudeImages(:,:,nn+1);
                        obj.setHolograms(intermediateField, intermediateHologramIndex);
                    end
                    % Propagate to object plane
                    propagationDistance = -obj.ZCoordinates(hologramCount);
                    propagatedWavelength = obj.Wavelengths(hologramCount);

                    obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                    obj.setWavelengths(propagatedWavelength, intermediateHologramIndex);

                    intermediateReconstruction = obj.propagate(intermediateHologramIndex);

                    % Rescale phase to 1st wavelength
                    phase = angle(intermediateReconstruction)*lambda(hologramCount)/lambda(1);
                    intermediateReconstruction = abs(intermediateReconstruction).*exp(1i.*phase);
                    obj.setHolograms(intermediateReconstruction, intermediateHologramIndex);

                    % Propagate to 1st hologram plane
                    propagationDistance = -obj.ZCoordinates(1);
                    propagatedWavelength = obj.Wavelengths(1);

                    obj.setZCoordinates(propagationDistance, intermediateHologramIndex);
                    obj.setWavelengths(propagatedWavelength, intermediateHologramIndex);

                    intermediateField = obj.propagate(intermediateHologramIndex);

                    % Actualize optical field with 1st amplitude
                    intermediateField = intermediateField./abs(intermediateField).*amplitudeImages(:,:,1);
                    obj.setHolograms(intermediateField, intermediateHologramIndex);
                end
            end

            % Backpropagate reconstructed optical field to object plane
            propagationDistance = -obj.ZCoordinates(1);
            obj.setZCoordinates(propagationDistance, intermediateHologramIndex);

            reconstruction = obj.propagate(intermediateHologramIndex);

            obj.ZCoordinates(intermediateHologramIndex) = [];
            obj.Holograms(:,:,intermediateHologramIndex) = [];

            if(sigma > 0)
                obj.Holograms = hologramsUnblurred;
            end
        end

        function reconstruction = GaborAveraging(obj)
            % Gabor averaging in-line holography reconstruction
            reconstruction = zeros(size(obj.Holograms(:,:,1)));
            hologramCount = size(obj.Holograms, 3);
            obj.ZCoordinates = -obj.ZCoordinates;

            for nn = 1:hologramCount
                % Backpropagate each hologram to object plane
                intermediateReconstruction = obj.propagate(nn);
                % Add the propagation result to the R_GA (and rescale phase to
                % match the 1st wavelength)
                reconstruction = reconstruction + sqrt(abs(intermediateReconstruction)) .* exp(1i.*angle(intermediateReconstruction) * obj.Wavelengths(nn)/obj.Wavelengths(1));
            end

            % Divide by the number of holograms
            reconstruction = reconstruction/hologramCount;

            obj.ZCoordinates = -obj.ZCoordinates;
        end
        % DarkFocus version 3
        function [bestFocusZ, focusCurves, DarkVolume] = DarkFocus(obj, hologramIndex, range, inspectedROI, darkVolumeROI)
            background = imgaussfilt(obj.Holograms(:,:,hologramIndex), 30);
            darkHologram = obj.Holograms(:,:,hologramIndex) - background;

            YIndexes = (inspectedROI(2):inspectedROI(2)+inspectedROI(4) - 1)+1;
            XIndexes = (inspectedROI(1):inspectedROI(1)+inspectedROI(3) - 1)+1;

            darkVolumeROI(2) = darkVolumeROI(2) - inspectedROI(2) + 1;
            darkVolumeROI(1) = darkVolumeROI(1) - inspectedROI(1) + 1;

            darkVolumeYIndexes = darkVolumeROI(2):(darkVolumeROI(2)+darkVolumeROI(4));
            darkVolumeXIndexes = darkVolumeROI(1):(darkVolumeROI(1)+darkVolumeROI(3));
            darkHologram = darkHologram(YIndexes,XIndexes);

            darkModel = QCI_Model("darkHologram",darkHologram,obj.Wavelengths(hologramIndex),obj.ZCoordinates(hologramIndex),obj.CameraPixelSize);
            DarkFocus = zeros(1,length(range));
            DarkVolume = zeros(length(darkVolumeYIndexes), length(darkVolumeXIndexes));

            for rangeIndex = 1:length(range)
                darkModel.setZCoordinates(range(rangeIndex),1);
                Obj1 = darkModel.propagate(1);
                % DarkFocus
                amp = abs(Obj1);
                DarkVolume(:,:,rangeIndex) = amp(darkVolumeYIndexes, darkVolumeXIndexes); %angle(Obj1);%
                [gx, gy] = gradient(DarkVolume(:,:,rangeIndex));
                grad = gx.^2+gy.^2;
                DarkFocus(rangeIndex) = var(grad(:));
            end
            
            % normalization 0-1
            DarkFocus = DarkFocus - min(DarkFocus); DarkFocus = DarkFocus/max(DarkFocus);
            focusCurves(1,:) = DarkFocus;
            [~,loc] = max(focusCurves);
            bestFocusZ = range(loc);

        end
    end
end

