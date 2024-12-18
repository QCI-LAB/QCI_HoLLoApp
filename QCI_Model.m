classdef QCI_Model < handle
    %QCI_MODEL is a bundle used for hologram creation and manipulation
    %   Detailed explanation goes here
    
    properties (Access = private)
        %Initial Data
        Names                   (1,:) string = []
        HologramsRaw            (:,:,:) double                              % Hologram image (matrix)
        Wavelengths             (1,:) double {mustBeNonnegative} = []       % Wavelength of light used [micro meters]
        ZCoordinates            (1,:) double {mustBeNonnegative} = []       % Distance from CCD to sample [micro meters] on the optical axis
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
        HologramsPostPropagation(:,:,:)                          = []
        
        CorrectedHolograms      (:,:,:)                          = []
        CorrectedReference      (:,:,:)                          = []
        CorrectedZCoordinates   (1,:) double {mustBeNonnegative} = []
    end
    
    methods
        function obj = QCI_Model(names, images, wavelengths, distances, pixelSize)
            if nargin == 0
                return
            end
            %QCI_MODEL Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Names = strings(1,length(names));
            obj.Names = names;
            for i = 1:length(wavelengths)
                obj.HologramsRaw(:,:,i) = double(cell2mat(images{i}));
            end
            obj.Wavelengths = wavelengths;
            obj.ZCoordinates = distances;
            obj.CameraPixelSize = pixelSize;
            [ySize, xSize] = size(obj.HologramsRaw(:,:,1));
            
            % Deriving data
            obj.YCoordinates = (0:(ySize-1))*obj.CameraPixelSize;
            obj.XCoordinates = (0:(xSize-1))*obj.CameraPixelSize;
            inverseStepY = 1/ySize/pixelSize;
            inverseStepX = 1/xSize/pixelSize;
            YInverseCoordinates =  -ySize/2*inverseStepY:inverseStepY:(ySize/2 - 1)*inverseStepY;
            XInverseCoordinates =  -xSize/2*inverseStepX:inverseStepX:(xSize/2 - 1)*inverseStepX;
            obj.WaveNumbers = 2*pi./obj.Wavelengths;
            for i = 1:length(wavelengths)
                obj.HologramsFFT(:,:,i) = fft2(obj.HologramsRaw(:,:,i));
                obj.FreeSpaceImpulseMatrixes(:,:,i) = fftshift(obj.WaveNumbers(i) *...
                    sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                    (ones(ySize,1)*(XInverseCoordinates.^2) + (YInverseCoordinates'.^2)*ones(1,xSize))));
            end
        end
        
        function names = getNames(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            names = obj.Names;
            return
        end

        function image = getImage(obj, hologramIndex)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            image = uint8(obj.HologramsRaw(:,:, hologramIndex));
            return
        end

        function image = getHologramPostPropagation(obj, hologramIndex)
            if(isempty(obj.HologramsPostPropagation))
                image = zeros(size(obj.HologramsRaw(:,:,1)));
                return
            end
            image = obj.HologramsPostPropagation(:,:,hologramIndex);
        end

        function image = getCorrectedHologram(obj, hologramIndex)
            if(isempty(obj.CorrectedHolograms))
                image = zeros(size(obj.HologramsRaw(:,:,1)));
                return
            end
            image = obj.CorrectedHolograms(:,:,hologramIndex);
        end

        function image = getCorrectedReference(obj, hologramIndex)
            if(isempty(obj.CorrectedReference))
                image = zeros(size(obj.HologramsPostPropagation(:,:,1)));
                return
            end
            
            image = obj.CorrectedReference(:,:,hologramIndex);
        end

        function Wavelengths = getWavelengths(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Wavelengths = obj.Wavelengths;
            return
        end

        function setWavelengths(obj, Wavelengths)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            for i = 1:length(Wavelengths)
                if(Wavelengths(i) ~= obj.Wavelengths(i))
                    obj.FreeSpaceImpulseMatrixes(:,:,i) = fftshift(obj.WaveNumbers(i) *...
                        sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                        (ones(ySize,1)*(XInverseCoordinates.^2) + (YInverseCoordinates'.^2)*ones(1,xSize))));
                end
            end
            obj.Wavelengths = Wavelengths;
            return
        end

        function ZCoordinates = getZCoordinates(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            ZCoordinates = obj.ZCoordinates;
            return
        end

        function setZCoordinates(obj, ZCoordinates)
            % Setter fro ZCoordinates
            obj.ZCoordinates = ZCoordinates;
            return
        end
        
        % Operations on holograms
        function hologramsPostPropagation = propagate(obj, z, hologramIndex)
            % Propagates the image with Angular Scaling method
            kernelExponent = obj.FreeSpaceImpulseMatrixes(:,:,hologramIndex)*z;
            kernelExponent = kernelExponent - kernelExponent(1,1);
            
            if(z<0)
                kernel = exp(-1i*kernelExponent);
                propagetedFFT = kernel.*obj.HologramsFFT(:,:, hologramIndex);
                obj.HologramsPostPropagation(:,:, hologramIndex) = conj(ifft2(propagetedFFT));
                hologramsPostPropagation = obj.HologramsPostPropagation(:,:, hologramIndex);
            else
                kernel = exp(1i*kernelExponent);
                propagetedFFT = kernel.*obj.HologramsFFT(:,:, hologramIndex);
                obj.HologramsPostPropagation(:,:, hologramIndex) = ifft2(propagetedFFT);
                hologramsPostPropagation = obj.HologramsPostPropagation(:,:, hologramIndex);
            end
        end

        function [fixedReferenceImages, fixedHolograms, tforms] = automaticalShiftAndScalingCorrection(obj, referenceImages, holograms)
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
            if (nargin == 0)
                if(~isempty(obj.HologramsRaw) && ~isempty(obj.HologramsPostPropagation))
                    referenceImages = obj.HologramsPostPropagation;
                    holograms = obj.HologramsRaw;
                else
                    error("Module has no holograms to work with.");
                end
            end

            imageNumber = size(referenceImages, 3);
            fixedReferenceImages = zeros(size(referenceImages), 'like', referenceImages);
            fixedHolograms = zeros(size(holograms), 'like', holograms);
            tforms = cell(1, imageNumber-1);  % Store transformations
        
            % Select the last reference image as the alignment target
            referenceImageFixed = referenceImages(:, :, end);
            fixedReferenceImages(:, :, end) = referenceImageFixed;
            fixedHolograms(:, :, end) = holograms(:, :, end);
            
            % Detect features in the fixed reference image
            featureDetectionThreshold = 10000;
            ptsFixed = detectSURFFeatures(referenceImageFixed, 'MetricThreshold', featureDetectionThreshold);
            [featuresFixed, validPtsFixed] = extractFeatures(referenceImageFixed, ptsFixed);
        
            % --- Loop through all images except the last one ---
            for k = 1:imageNumber-1
                % Current distorted reference image and hologram
                distortedImage = referenceImages(:, :, k);
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
        
                % Optional Debugging (uncomment for visualization)
                % figure; showMatchedFeatures(referenceImageFixed, distortedImage, ...
                %                            matchedFixed(inlierIdx), matchedDistorted(inlierIdx));
                % title(sprintf('Image %d Matching Points (Inliers Only)', k));
            end
            obj.CorrectedHolograms = fixedHolograms;
            obj.CorrectedReference = fixedReferenceImages;
        end

        function Reconstruction = IGA(obj, Hologram,z,wavelength,dx,sigma,iter)
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
            
            [Sy,Sx,N] = size(Hologram);
            
            if sigma < inf % if not GA method
                if sigma > 0
                    % Low-pass filtering of input holograms
                    Hblurred = zeros(Sy,Sx,N);
                    for n = 1:N
                        Hblurred(:,:,n) = imgaussfilt(Hologram(:,:,n),sigma);
                    end
                else
                    Hblurred = Hologram; % No bluring for GS method
                end
                % Gerchberg-Saxton reconstruction
                if isscalar(length(unique(obj.Wavelengths))) 
                    R_GS = GS_multiHeight(Hblurred, z, obj.Wavelengths(1),dx,iter,N);
                else
                    R_GS = GS_multiWavelength(Hblurred, z, obj.Wavelengths(1),dx,iter,N);
                end
            end
            
            % Gabor averaging reconstruction
            if sigma > 0 % if not the GS method
                R_GA = GA(Hologram,z,wavelength,dx,N);
            end
            
            if sigma == 0 % GS method
                Reconstruction = R_GS;
            elseif sigma == inf % GA method
                Reconstruction = R_GA;
            else % IGA method
                % Combine GS and GA methods
                Rea = real(R_GS) + real(R_GA) - imgaussfilt(real(R_GA),sigma);
                Ima = imag(R_GS) + imag(R_GA) - imgaussfilt(imag(R_GA),sigma);
                Reconstruction = Rea + 1i.*Ima;
            end
        end

        function Reconstruction = GaborAveraging(obj, H, z)
        % Gabor averaging in-line holography reconstruction
            Reconstruction = 0;
            hologramNumber = size(obj.HologramsRaw, 3);
            for nn = 1:hologramNumber
                % Backpropagate each hologram to object plane
                U = AS_propagate_p(H(:,:,nn), -obj.ZCoordinates(nn), obj.Wavelengths(nn), obj.CameraPixelSize);
                % Add the propagation result to the R_GA (and rescale phase to
                % match the 1st wavelength)
                Reconstruction = Reconstruction + sqrt(abs(U)) .* exp(1i.*angle(U) * obj.Wavelengths(nn)/obj.Wavelengths(1));
            end
            % Divide by the number of holograms
            Reconstruction = Reconstruction/hologramNumber;
        end

    end
end

