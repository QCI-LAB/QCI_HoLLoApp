classdef QCI_Model < handle
    %QCI_MODEL is a bundle used for hologram creation and manipulation
    %   Detailed explanation goes here
    
    properties (Access = protected)
        %Initial Data
        Names                   (1,:) string = []
        Holograms               (:,:,:) double                   = []       % Hologram image (matrix)
        Wavelengths             (1,:) double {mustBeNonnegative} = []       % Wavelength of light used [micro meters]
        PropagationDistances    (1,:) double                     = []       % Distance from CCD to sample [micro meters] on the optical axis
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
        function obj = QCI_Model(names, holograms, wavelengths, propagationDistances, pixelSize)
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
            obj.PropagationDistances = propagationDistances;
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
            obj.PropagationDistances(end + 1) = zCoordinate;
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
            obj.PropagationDistances(hologramIndex) = [];
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
            [ySize, xSize] = size(obj.Holograms(:,:,1));
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
                for i = 1:length(indexes)
                    if(wavelengths(i) ~= obj.Wavelengths(indexes(i)))
                        obj.FreeSpaceImpulseMatrixes(:,:,indexes(i)) = fftshift(obj.WaveNumbers(indexes(i)) *...
                            sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(indexes(i))^2 *...
                            (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
                    end
                end

                obj.Wavelengths(indexes) = wavelengths;
            else
                error("Dimension of indexes array does not match the amount of given wavelengths.");
            end

            return
        end

        function propagationDistances = getPropagationDistances(obj, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if(indexes == 0)
                propagationDistances = obj.PropagationDistances;
            else
                propagationDistances = obj.PropagationDistances(indexes);
            end
        end

        function setPropagationDistances(obj, propagationDistances, indexes)
            % Setter for ZCoordinates
            if(size(propagationDistances) == size(obj.PropagationDistances))
                obj.PropagationDistances = propagationDistances;
            elseif(size(propagationDistances) == size(indexes))  
                obj.PropagationDistances(indexes) = propagationDistances;
            else
                error("Dimension of indexes array does not match the amount of given propagationDistances.");
            end
        end

        function cameraPixelSize = getCameraPixelSize(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cameraPixelSize = obj.CameraPixelSize;
        end

        function setCameraPixelSize(obj, newCameraPixelSize)
            % Setter for ZCoordinates
            if(newCameraPixelSize > 0)
                obj.CameraPixelSize = newCameraPixelSize;
            else
                error("Camera pixel size must be a non-negative value.");
            end
        end
        
        % Operations on holograms
        function hologramPostPropagation = propagate(obj, hologramIndex, willReplace)
            % Propagates the image with Angular Scaling method
            kernelExponent = obj.FreeSpaceImpulseMatrixes(:,:,hologramIndex)*obj.PropagationDistances(hologramIndex);
            kernelExponent = kernelExponent - kernelExponent(1,1);
            %hologramPostPropagation = zeros(size(obj.Holograms,1), size(obj.Holograms,2));

            if(obj.PropagationDistances(hologramIndex)<0)
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
                correctedHeights(tt) = round(obj.PropagationDistances(tt).*tform.T(1,1).^2);
                tforms{tt} = tform;
            end
            fixedHolograms(:,:,end+1) = obj.Holograms(:,:,end);
            fixedReferenceImages(:,:,end+1) = referenceImages(:,:,end);
            correctedHeights(end+1) = obj.PropagationDistances(end);
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
            correctedHeights(end) = obj.PropagationDistances(end);
            
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
                correctedHeights(k) = round(obj.PropagationDistances(k).*tform.T(1,1).^2);
        
                % Optional Debugging (uncomment for visualization)
                % figure; showMatchedFeatures(referenceImageFixed, distortedImage, ...
                %                            matchedFixed(inlierIdx), matchedDistorted(inlierIdx));
                % title(sprintf('Image %d Matching Points (Inliers Only)', k));
            end
        end

        function  fixedHolograms = automaticalShiftCorrection(obj, holograms)
            % --- Initialization ---
            if(isempty(holograms))
                error("No holograms to work with.");
            end

            hologramsFFT = fft2(holograms);

            usfactor = 10;

            for hologramIndex = 1:(size(holograms,3)-1)
                [~, Greg] = dftregistration(hologramsFFT(:,:,end), hologramsFFT(:,:,hologramIndex), usfactor)
                fixedHolograms(:,:,hologramIndex) = ifft2(Greg);
            end
            fixedHolograms(:,:,end) = holograms(:,:,end);
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

            obj.addHologram("intermediateHologram", amplitudeImages(:, :, 1), obj.Wavelengths(1), obj.PropagationDistances(1)); % ZCoordinate will be replaced
            intermediateHologramIndex = hologramCount + 1;

            if(isSingleWavelength)
                for tt = 1:iter
                    for nn = 1:(hologramCount-1)

                        % Propagate to nn+1 hologram plane
                        propagationDistance = obj.PropagationDistances(nn+1) - obj.PropagationDistances(nn);
                        obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                        intermediateField = obj.propagate(intermediateHologramIndex);

                        % Actualize optical field with nn+1 amplitude
                        intermediateField = intermediateField./abs(intermediateField) .* amplitudeImages(:,:,nn+1);
                        obj.setHolograms(intermediateField, intermediateHologramIndex);
                    end
                    % Propagate to 1st hologram plane
                    propagationDistance = obj.PropagationDistances(1) - obj.PropagationDistances(hologramCount);
                    obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                    intermediateField = obj.propagate(intermediateHologramIndex);
                    % Actualize optical field with 1st amplitude
                    intermediateField = intermediateField./abs(intermediateField) .* amplitudeImages(:,:,1);
                    obj.setHolograms(intermediateField, intermediateHologramIndex);
                end
            else
                for tt = 1:iter
                    for nn = 1:(hologramCount-1)
                        % Propagate to object plane
                        propagationDistance = -obj.PropagationDistances(nn);
                        obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                        intermediateReconstruction = obj.propagate(intermediateHologramIndex);

                        % Rescale phase to nn+1 wavelength
                        phase = angle(intermediateReconstruction) * obj.Wavelengths(nn)/obj.Wavelengths(nn+1);
                        intermediateReconstruction = abs(intermediateReconstruction) .* exp(1i.*phase);
                        obj.setHolograms(intermediateReconstruction, intermediateHologramIndex);

                        % Propagate to nn+1 hologram plane
                        propagationDistance = obj.PropagationDistances(nn+1);
                        obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                        obj.setWavelengths(obj.Wavelengths(nn+1), intermediateHologramIndex);
                        intermediateField = obj.propagate(intermediateHologramIndex);

                        % Actualize optical field with nn+1 amplitude
                        intermediateField = intermediateField ./ abs(intermediateField) .* amplitudeImages(:,:,nn+1);
                        obj.setHolograms(intermediateField, intermediateHologramIndex);
                    end
                    % Propagate to object plane
                    propagationDistance = -obj.PropagationDistances(hologramCount);
                    propagatedWavelength = obj.Wavelengths(hologramCount);

                    obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                    obj.setWavelengths(propagatedWavelength, intermediateHologramIndex);

                    intermediateReconstruction = obj.propagate(intermediateHologramIndex);

                    % Rescale phase to 1st wavelength
                    phase = angle(intermediateReconstruction)*obj.Wavelengths(hologramCount)/obj.Wavelengths(1);
                    intermediateReconstruction = abs(intermediateReconstruction).*exp(1i.*phase);
                    obj.setHolograms(intermediateReconstruction, intermediateHologramIndex);

                    % Propagate to 1st hologram plane
                    propagationDistance = -obj.PropagationDistances(1);
                    propagatedWavelength = obj.Wavelengths(1);

                    obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);
                    obj.setWavelengths(propagatedWavelength, intermediateHologramIndex);

                    intermediateField = obj.propagate(intermediateHologramIndex);

                    % Actualize optical field with 1st amplitude
                    intermediateField = intermediateField./abs(intermediateField).*amplitudeImages(:,:,1);
                    obj.setHolograms(intermediateField, intermediateHologramIndex);
                end
            end

            % Backpropagate reconstructed optical field to object plane
            propagationDistance = -obj.PropagationDistances(1);
            obj.setPropagationDistances(propagationDistance, intermediateHologramIndex);

            reconstruction = obj.propagate(intermediateHologramIndex);

            obj.PropagationDistances(intermediateHologramIndex) = [];
            obj.Holograms(:,:,intermediateHologramIndex) = [];

            if(sigma > 0)
                obj.Holograms = hologramsUnblurred;
            end
        end

        function reconstruction = GaborAveraging(obj)
            % Gabor averaging in-line holography reconstruction
            reconstruction = zeros(size(obj.Holograms(:,:,1)));
            hologramCount = size(obj.Holograms, 3);
            obj.PropagationDistances = -obj.PropagationDistances;

            for nn = 1:hologramCount
                % Backpropagate each hologram to object plane
                intermediateReconstruction = obj.propagate(nn);
                % Add the propagation result to the R_GA (and rescale phase to
                % match the 1st wavelength)
                reconstruction = reconstruction + sqrt(abs(intermediateReconstruction)) .* exp(1i.*angle(intermediateReconstruction) * obj.Wavelengths(nn)/obj.Wavelengths(1));
            end

            % Divide by the number of holograms
            reconstruction = reconstruction/hologramCount;

            obj.PropagationDistances = -obj.PropagationDistances;
        end
        % DarkFocus version 3
        function [bestFocusZ, focusCurves, DarkVolume] = DarkFocus(obj, hologramIndex, range, inspectedROI, darkVolumeROI)
            background = imgaussfilt(obj.Holograms(:,:,hologramIndex), 30);
            darkHologram = obj.Holograms(:,:,hologramIndex) - background;

            YIndexes = (inspectedROI(2):inspectedROI(2)+inspectedROI(4) - 1) + 1;
            XIndexes = (inspectedROI(1):inspectedROI(1)+inspectedROI(3) - 1) + 1;
            if(darkVolumeROI == 0)
                darkVolumeROI = inspectedROI;
            else
                darkVolumeROI(2) = darkVolumeROI(2) - inspectedROI(2) + 1;
                darkVolumeROI(1) = darkVolumeROI(1) - inspectedROI(1) + 1;
            end
            darkVolumeYIndexes = (darkVolumeROI(2):(darkVolumeROI(2)+darkVolumeROI(4)) - 1) + 1;
            darkVolumeXIndexes = (darkVolumeROI(1):(darkVolumeROI(1)+darkVolumeROI(3)) - 1) + 1;
            darkHologram = darkHologram(YIndexes,XIndexes);

            darkModel = QCI_Model("darkHologram",darkHologram,obj.Wavelengths(hologramIndex),obj.PropagationDistances(hologramIndex),obj.CameraPixelSize);
            DarkFocus = zeros(1,length(range));
            DarkVolume = zeros(length(darkVolumeYIndexes), length(darkVolumeXIndexes));

            for rangeIndex = 1:length(range)
                darkModel.setPropagationDistances(range(rangeIndex),1);
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

        function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
        % function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
        % Efficient subpixel image registration by crosscorrelation. This code
        % gives the same precision as the FFT upsampled cross correlation in a
        % small fraction of the computation time and with reduced memory 
        % requirements. It obtains an initial estimate of the crosscorrelation peak
        % by an FFT and then refines the shift estimation by upsampling the DFT
        % only in a small neighborhood of that estimate by means of a 
        % matrix-multiply DFT. With this procedure all the image points are used to
        % compute the upsampled crosscorrelation.
        % Manuel Guizar - Dec 13, 2007
        %
        % Rewrote all code not authored by either Manuel Guizar or Jim Fienup
        % Manuel Guizar - May 13, 2016
        %
        % Citation for this algorithm:
        % Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
        % "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
        % 156-158 (2008).
        %
        % Inputs
        % buf1ft    Fourier transform of reference image, 
        %           DC in (1,1)   [DO NOT FFTSHIFT]
        % buf2ft    Fourier transform of image to register, 
        %           DC in (1,1) [DO NOT FFTSHIFT]
        % usfac     Upsampling factor (integer). Images will be registered to 
        %           within 1/usfac of a pixel. For example usfac = 20 means the
        %           images will be registered within 1/20 of a pixel. (default = 1)
        %
        % Outputs
        % output =  [error,diffphase,net_row_shift,net_col_shift]
        % error     Translation invariant normalized RMS error between f and g
        % diffphase     Global phase difference between the two images (should be
        %               zero if images are non-negative).
        % net_row_shift net_col_shift   Pixel shifts between images
        % Greg      (Optional) Fourier transform of registered version of buf2ft,
        %           the global phase difference is compensated for.
        %
        %
        % Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
        % All rights reserved.
        % 
        % Redistribution and use in source and binary forms, with or without
        % modification, are permitted provided that the following conditions are
        % met:
        % 
        %     * Redistributions of source code must retain the above copyright
        %       notice, this list of conditions and the following disclaimer.
        %     * Redistributions in binary form must reproduce the above copyright
        %       notice, this list of conditions and the following disclaimer in
        %       the documentation and/or other materials provided with the distribution
        %     * Neither the name of the University of Rochester nor the names
        %       of its contributors may be used to endorse or promote products derived
        %       from this software without specific prior written permission.
        % 
        % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
        % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
        % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
        % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
        % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
        % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
        % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
        % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
        % POSSIBILITY OF SUCH DAMAGE.
            
            if ~exist('usfac','var')
                usfac = 1;
            end
            [nr,nc]=size(buf2ft);
            Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            if usfac == 0
                % Simple computation of error and phase difference without registration
                CCmax = sum(buf1ft(:).*conj(buf2ft(:)));
                row_shift = 0;
                col_shift = 0;
            elseif usfac == 1
                % Single pixel registration
                CC = ifft2(buf1ft.*conj(buf2ft));
                CCabs = abs(CC);
                [row_shift, col_shift] = find(CCabs == max(CCabs(:)));
                CCmax = CC(row_shift,col_shift)*nr*nc;
                % Now change shifts so that they represent relative shifts and not indices
                row_shift = Nr(row_shift);
                col_shift = Nc(col_shift);
            elseif usfac > 1
                % Start with usfac == 2
                CC = ifft2(FTpad(buf1ft.*conj(buf2ft),[2*nr,2*nc]));
                CCabs = abs(CC);
                [row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');
                CCmax = CC(row_shift,col_shift)*nr*nc;
                % Now change shifts so that they represent relative shifts and not indices
                Nr2 = ifftshift(-fix(nr):ceil(nr)-1);
                Nc2 = ifftshift(-fix(nc):ceil(nc)-1);
                row_shift = Nr2(row_shift)/2;
                col_shift = Nc2(col_shift)/2;
                % If upsampling > 2, then refine estimate with matrix multiply DFT
                if usfac > 2
                    %%% DFT computation %%%
                    % Initial shift estimate in upsampled grid
                    row_shift = round(row_shift*usfac)/usfac; 
                    col_shift = round(col_shift*usfac)/usfac;     
                    dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
                    % Matrix multiply DFT around the current shift estimate
                    CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
                        dftshift-row_shift*usfac,dftshift-col_shift*usfac));
                    % Locate maximum and map back to original pixel grid 
                    CCabs = abs(CC);
                    [rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');
                    CCmax = CC(rloc,cloc);
                    rloc = rloc - dftshift - 1;
                    cloc = cloc - dftshift - 1;
                    row_shift = row_shift + rloc/usfac;
                    col_shift = col_shift + cloc/usfac;    
                end
                % If its only one row or column the shift along that dimension has no
                % effect. Set to zero.
                if nr == 1
                    row_shift = 0;
                end
                if nc == 1
                    col_shift = 0;
                end
                
            end  
            rg00 = sum(abs(buf1ft(:)).^2);
            rf00 = sum(abs(buf2ft(:)).^2);
            error = 1.0 - abs(CCmax).^2/(rg00*rf00);
            error = sqrt(abs(error));
            diffphase = angle(CCmax);
            output=[error,diffphase,row_shift,col_shift];
            % Compute registered version of buf2ft
            if (nargout > 1)&&(usfac > 0)
                [Nc,Nr] = meshgrid(Nc,Nr);
                Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
                Greg = Greg*exp(1i*diffphase);
            elseif (nargout > 1)&&(usfac == 0)
                Greg = buf2ft*exp(1i*diffphase);
            end
        end
        function out=dftups(in,nor,noc,usfac,roff,coff)
        % function out=dftups(in,nor,noc,usfac,roff,coff);
        % Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
        % a small region.
        % usfac         Upsampling factor (default usfac = 1)
        % [nor,noc]     Number of pixels in the output upsampled DFT, in
        %               units of upsampled pixels (default = size(in))
        % roff, coff    Row and column offsets, allow to shift the output array to
        %               a region of interest on the DFT (default = 0)
        % Recieves DC in upper left corner, image center must be in (1,1) 
        % Manuel Guizar - Dec 13, 2007
        % Modified from dftus, by J.R. Fienup 7/31/06
        % This code is intended to provide the same result as if the following
        % operations were performed
        %   - Embed the array "in" in an array that is usfac times larger in each
        %     dimension. ifftshift to bring the center of the image to (1,1).
        %   - Take the FFT of the larger array
        %   - Extract an [nor, noc] region of the result. Starting with the 
        %     [roff+1 coff+1] element.
        % It achieves this result by computing the DFT in the output array without
        % the need to zeropad. Much faster and memory efficient than the
        % zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]
            [nr,nc]=size(in);
            % Set defaults
            if exist('roff', 'var')~=1, roff=0;  end
            if exist('coff', 'var')~=1, coff=0;  end
            if exist('usfac','var')~=1, usfac=1; end
            if exist('noc',  'var')~=1, noc=nc;  end
            if exist('nor',  'var')~=1, nor=nr;  end
            % Compute kernels and obtain DFT by matrix products
            kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
            kernr=exp((-1i*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
            out=kernr*in*kernc;
        end
        
        function [ imFTout ] = FTpad(imFT,outsize)
            % imFTout = FTpad(imFT,outsize)
            % Pads or crops the Fourier transform to the desired ouput size. Taking 
            % care that the zero frequency is put in the correct place for the output
            % for subsequent FT or IFT. Can be used for Fourier transform based
            % interpolation, i.e. dirichlet kernel interpolation. 
            %
            %   Inputs
            % imFT      - Input complex array with DC in [1,1]
            % outsize   - Output size of array [ny nx] 
            %
            %   Outputs
            % imout   - Output complex image with DC in [1,1]
            % Manuel Guizar - 2014.06.02
            if ~ismatrix(imFT)
                error('Maximum number of array dimensions is 2')
            end
            Nout = outsize;
            Nin = size(imFT);
            imFT = fftshift(imFT);
            center = floor(size(imFT)/2)+1;
            imFTout = zeros(outsize);
            centerout = floor(size(imFTout)/2)+1;
            % imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
            %     = imFT;
            cenout_cen = centerout - center;
            imFTout( ...
                max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)), ...
                max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2)) ...
                ) ...
                = imFT( ...
                max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1), ...
                Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)) ...
                );
            imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)/(Nin(1)*Nin(2));
        end
    end
end

