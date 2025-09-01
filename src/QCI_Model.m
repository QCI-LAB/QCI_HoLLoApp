classdef QCI_Model < handle
    %QCI_MODEL is a bundle used for hologram creation and manipulation
    %   Detailed explanation goes here
    
    properties (Access = protected)
        %Initial Data
        Names                   (1,:) string                     = []
        Holograms               (:,:,:) double                   = []       % Hologram image (matrix)
        Wavelengths             (1,:) double {mustBeNonnegative} = []       % Wavelength of light used [micro meters]
        PropagationDistances    (1,:) double                     = []       % Distance from CCD to sample [micro meters] on the optical axis
        CameraPixelSize         (1,1) double                     = 1        % Size of a pixel on the CCD camera used to obtain the hologram
        MediumRefractiveIndex   (1,1) double                     = 1        % Refractive index of the medium in which the hologram was taken (assume air)

        % Derrived data
        YCoordinates            (1,:) double {mustBeNonnegative} = []       % Y hologram coordinates in physical units [micro meters]. Y axis points down from the left upper corner.
        XCoordinates            (1,:) double {mustBeNonnegative} = []       % X hologram coordinates in physical units [micro meters]. X axis points right from the left upper corner.
        YInverseCoordinates     (1,:) double                     = []       % Y hologram coordinates in inverse space;
        XInverseCoordinates     (1,:) double                     = []       % X hologram coordinates in inverse space;
        WaveNumbers             (1,:) double {mustBeNonnegative} = []       % Wavenumber in [1/micrometer = 1e6 m]
        PropagationKernels      (:,:,:) double                   = []
        HologramsFFT            (:,:,:)                          = []
    end
    
    methods
        function obj = QCI_Model(names, holograms, wavelengths, propagationDistances, pixelSize)
            % QCI_Model Constructor for the QCI_Model class  
            %  
            %   obj = QCI_Model(names, holograms, wavelengths, propagationDistances, pixelSize)  
            %   initializes an instance with hologram data and computes derived parameters.  
            %  
            %   Inputs:  
            %       - names: (1×N string) Hologram names.  
            %       - holograms: (M×N×K double) Hologram images.  
            %       - wavelengths: (1×K double) Light wavelengths [μm].  
            %       - propagationDistances: (1×K double) Axial distances [μm].  
            %       - pixelSize: (scalar double) Camera pixel size [μm].  
            %  
            %   Computes spatial and inverse coordinates, wavenumbers, FFT of holograms,  
            %   and free-space impulse matrices for hologram reconstruction.  
            %   Initializes an empty object if no inputs are provided.

            if nargin == 0
                return
            end
            
            obj.Names = names;
            hologramCount = size(holograms,3);
            for i = 1:hologramCount
                obj.Holograms(:,:,i) = holograms(:,:,i);
            end

            obj.Wavelengths = wavelengths;
            obj.PropagationDistances = propagationDistances;
            obj.CameraPixelSize = pixelSize;
            [ySize, xSize] = size(obj.Holograms(:,:,1));
            
            % Deriving data
            obj.YCoordinates = (0:(ySize-1))*obj.CameraPixelSize;
            obj.XCoordinates = (0:(xSize-1))*obj.CameraPixelSize;
            inverseStepY = 1/ySize/obj.CameraPixelSize;
            inverseStepX = 1/xSize/obj.CameraPixelSize;
            obj.YInverseCoordinates =  -ySize/2*inverseStepY:inverseStepY:(ySize/2 - 1)*inverseStepY;
            obj.XInverseCoordinates =  -xSize/2*inverseStepX:inverseStepX:(xSize/2 - 1)*inverseStepX;
            obj.WaveNumbers = 2*pi./obj.Wavelengths;

            obj.PropagationKernels = zeros(ySize, xSize, hologramCount);
            
            for i = 1:hologramCount
                obj.HologramsFFT(:,:,i) = fft2(obj.Holograms(:,:,i));
                phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                    sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                    (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                phi = phi - phi(1,1);
                obj.PropagationKernels(:,:,i) = exp(1i*phi);
            end
        end
        
        function names = getNames(obj, indexes)
            % getNames Returns hologram names for given indexes  
            %  
            %   names = getNames(obj, indexes) returns the names of the holograms  
            %   corresponding to the specified indexes.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - indexes: (1×N double) Indices of the requested holograms.  
            %         Use 0 to return all names.  
            %  
            %   Output:  
            %       - names: (1×N string) Selected hologram names.  
            if(indexes == 0)
                names = obj.Names;
            else
                names = obj.Names(indexes);
            end
     
            return
        end

        function setNames(obj, names, indexes)
            % setNames Updates the Names property for the given indexes
            %
            %   setNames(obj, names, indexes) updates the Names property of the object.
            %   If indexes is 0, the entire Names property is replaced with the provided names.
            %   Otherwise, only the entries at the specified indexes are updated.
            %
            %   Inputs:
            %       - obj: Instance of the class.
            %       - names: (e.g., cell array or string array) The new names to be set.
            %       - indexes: (1×N double) Indices where the names should be updated.
            %         Use 0 to replace all names.
            %
            %   Output:
            %       None (the method updates the object's Names property directly).
            
            if(indexes == 0)
                obj.Names = names;
            else
                obj.Names(indexes) = names;
            end
     
            return
        end

        function holograms = getHolograms(obj, indexes)
            % getHolograms Returns hologram images for given indexes  
            %  
            %   holograms = getHolograms(obj, indexes) returns the hologram images  
            %   corresponding to the specified indexes.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - indexes: (1×N double) Indices of the requested holograms.  
            %         Use 0 to return all holograms.  
            %  
            %   Output:  
            %       - holograms: (M×N×K double) Selected hologram images.

            if(indexes == 0)
                holograms = obj.Holograms;
            else
                holograms = obj.Holograms(:,:, indexes);
            end
        end
        

        function setHolograms(obj, holograms, indexes)
            % setHolograms Updates hologram images and their FFT  
            %  
            %   setHolograms(obj, holograms, indexes) updates the hologram images  
            %   in the object. The function checks if the size of the provided  
            %   holograms matches the expected dimensions and updates the hologram  
            %   data and its FFT.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - holograms: (M×N×K double) New hologram images.  
            %       - indexes: (1×N double) Indices for which holograms should be updated.  
            %  
            %   Throws an error if the number of given holograms does not match the  
            %   number of specified indexes.  

            if(size(holograms, 3) == size(obj.Holograms,3))
                obj.Holograms = double(holograms);
                obj.HologramsFFT = fft2(holograms);
            elseif(size(holograms, 3) == length(indexes))
                obj.Holograms(:,:,indexes) = double(holograms);
                obj.HologramsFFT(:,:,indexes) = fft2(holograms);
            else
                error("Dimension of indexes array does not match the amount of given holograms.");
            end
        end

        function addHologram(obj, name, hologram, wavelength, zCoordinate)
            % addHologram Adds a new hologram to the model  
            %  
            %   addHologram(obj, name, hologram, wavelength, zCoordinate) adds a new hologram  
            %   to the model with the specified name, hologram data, wavelength, and propagation  
            %   distance. It also computes the hologram's FFT and updates the corresponding  
            %   free-space impulse matrix.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - name: (string) Name of the new hologram.  
            %       - hologram: (M×N double) Hologram image data.  
            %       - wavelength: (double) Wavelength of the light used [μm].  
            %       - zCoordinate: (double) Propagation distance [μm].

            obj.Names(end + 1) = name;
            obj.Holograms(:, :, end + 1) = hologram;
            obj.Wavelengths(end + 1) = wavelength;
            obj.PropagationDistances(end + 1) = zCoordinate;
            obj.WaveNumbers(end + 1) = 2*pi./obj.Wavelengths(end);
            
            obj.HologramsFFT(:,:,end + 1) = fft2(obj.Holograms(:,:,end));
            [ySize, xSize] = size(obj.Holograms(:, :, 1));

            phi = real(fftshift(obj.WaveNumbers(end + 1) *abs(propagationDistances(end + 1))*...
                    sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(end + 1)^2 *...
                    (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
            phi = phi - phi(1,1);
            obj.PropagationKernels(:,:,end + 1) = exp(1i*phi);
        end

        function removeHologram(obj, hologramIndex)
            % removeHologram Removes a hologram from the model  
            %  
            %   removeHologram(obj, hologramIndex) removes the hologram at the specified index  
            %   and updates the model's data accordingly. If all holograms are removed,  
            %   it resets the coordinates and other related properties.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - hologramIndex: (scalar double) Index of the hologram to be removed. 
            obj.Holograms(:, :, hologramIndex) = [];
            obj.Wavelengths(hologramIndex) = [];
            obj.PropagationDistances(hologramIndex) = [];
            obj.WaveNumbers(hologramIndex) = [];
            
            obj.HologramsFFT(:, :, hologramIndex) = [];
            obj.PropagationKernels(:, :, hologramIndex) = [];

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
            % getWavelengths Returns wavelengths for given indexes  
            %  
            %   wavelengths = getWavelengths(obj, indexes) returns the wavelengths  
            %   corresponding to the specified indexes.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - indexes: (1×N double) Indices of the requested wavelengths.  
            %         Use 0 to return all wavelengths.  
            %  
            %   Output:  
            %       - wavelengths: (1×N double) Selected wavelengths.

            if(indexes == 0)
                wavelengths = obj.Wavelengths;
            else
                wavelengths = obj.Wavelengths(indexes);
            end
        end

        function setWavelengths(obj, wavelengths, indexes)
            % setWavelengths Updates wavelengths and related parameters  
            %  
            %   setWavelengths(obj, wavelengths, indexes) updates the wavelengths and  
            %   computes the corresponding wavenumbers and free-space impulse matrices.  
            %   The function checks if the provided wavelengths match the existing ones  
            %   and updates the relevant data accordingly.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - wavelengths: (1×N double) New wavelengths to be set.  
            %       - indexes: (1×N double) Indices of the wavelengths to be updated.  
            %  
            %   Throws an error if the number of given wavelengths does not match the  
            %   number of specified indexes. 

            [ySize, xSize] = size(obj.Holograms(:,:,1));
            if(length(wavelengths) == length(obj.Wavelengths))
                for i = 1:length(wavelengths)
                    if(wavelengths(i) ~= obj.Wavelengths(i))
                        obj.Wavelengths(i) = wavelengths(i);
                        obj.WaveNumbers(i) = 2*pi/wavelengths(i);
                        phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                                    sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                                    (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                        phi = phi - phi(1,1);
                        obj.PropagationKernels(:,:,i) = exp(1i*phi);
                    end
                end
                obj.Wavelengths = wavelengths;
            elseif(length(wavelengths) == length(indexes))
                obj.Wavelengths(indexes) = wavelengths;
                obj.WaveNumbers(indexes) = 2*pi/wavelengths;
                
                for i = 1:length(indexes)       
                    phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                                sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                                (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                    phi = phi - phi(1,1);
                    obj.PropagationKernels(:,:,i) = exp(1i*phi);
                end
            else
                error("Dimension of indexes array does not match the amount of given wavelengths.");
            end

            return
        end

        function propagationDistances = getPropagationDistances(obj, indexes)
            % getPropagationDistances Returns propagation distances for given indexes  
            %  
            %   propagationDistances = getPropagationDistances(obj, indexes) returns the  
            %   propagation distances corresponding to the specified indexes.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - indexes: (1×N double) Indices of the requested propagation distances.  
            %         Use 0 to return all propagation distances.  
            %  
            %   Output:  
            %       - propagationDistances: (1×N double) Selected propagation distances.

            if(indexes == 0)
                propagationDistances = obj.PropagationDistances;
            else
                propagationDistances = obj.PropagationDistances(indexes);
            end
        end

        function setPropagationDistances(obj, propagationDistances, indexes)
            % setPropagationDistances Updates propagation distances  
            %  
            %   setPropagationDistances(obj, propagationDistances, indexes) updates the  
            %   propagation distances for the specified indexes. If the number of provided  
            %   distances matches the existing ones, they are updated.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - propagationDistances: (1×N double) New propagation distances.  
            %       - indexes: (1×N double) Indices of the distances to be updated.  
            %  
            %   Throws an error if the number of given propagation distances does not match  
            %   the number of specified indexes.
            [ySize, xSize] = size(obj.Holograms(:,:,1));
            if(length(propagationDistances) == length(obj.PropagationDistances))
                for i =1:length(propagationDistances)
                    if(propagationDistances(i)~=obj.PropagationDistances(i))
                        obj.PropagationDistances(i) = propagationDistances(i);
                        phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                            sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                            (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                        phi = phi - phi(1,1);
                        obj.PropagationKernels(:,:,i) = exp(1i*phi);
                    end
                end
            elseif(length(propagationDistances) == length(indexes))
                for i = indexes
                    if(propagationDistances(i)~=obj.PropagationDistances(i))
                        obj.PropagationDistances(i) = propagationDistances(i);
                        phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                            sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                            (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                        phi = phi - phi(1,1);
                        obj.PropagationKernels(:,:,i) = exp(1i*phi);
                    end
                end
            else
                error("Dimension of indexes array does not match the amount of given propagationDistances.");
            end
        end

        function cameraPixelSize = getCameraPixelSize(obj)
            % getCameraPixelSize Returns the camera pixel size  
            %  
            %   cameraPixelSize = getCameraPixelSize(obj) returns the size of the camera  
            %   pixels used in the hologram acquisition.  
            %  
            %   Output:  
            %       - cameraPixelSize: (scalar double) Size of the camera pixel [μm].

            cameraPixelSize = obj.CameraPixelSize;
        end

        function setCameraPixelSize(obj, newCameraPixelSize)
            % setCameraPixelSize Sets a new camera pixel size  
            %  
            %   setCameraPixelSize(obj, newCameraPixelSize) sets the camera pixel size  
            %   to the specified value. The new value must be positive.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - newCameraPixelSize: (scalar double) New camera pixel size [μm].  
            %  
            %   Throws an error if the new camera pixel size is non-positive.
            if(newCameraPixelSize > 0)
                [ySize, xSize] = size(obj.Holograms(:,:,1));
                obj.CameraPixelSize = newCameraPixelSize;
                % Deriving data
                obj.YCoordinates = (0:(ySize-1))*obj.CameraPixelSize;
                obj.XCoordinates = (0:(xSize-1))*obj.CameraPixelSize;
                inverseStepY = 1/ySize/obj.CameraPixelSize;
                inverseStepX = 1/xSize/obj.CameraPixelSize;
                obj.YInverseCoordinates =  -ySize/2*inverseStepY:inverseStepY:(ySize/2 - 1)*inverseStepY;
                obj.XInverseCoordinates =  -xSize/2*inverseStepX:inverseStepX:(xSize/2 - 1)*inverseStepX;
                obj.WaveNumbers = 2*pi./obj.Wavelengths;
                
                for i = 1:length(obj.Wavelengths)
                    phi = real(fftshift(obj.WaveNumbers(i) *abs(propagationDistances(i))*...
                                sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(i)^2 *...
                                (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize)))));
                    phi = phi - phi(1,1);
                    obj.PropagationKernels(:,:,i) = exp(1i*phi);
                end
            else
                error("Camera pixel size must be a non-negative value.");
            end
        end
        
        % Operations on holograms
        function hologramPostPropagation = propagate(obj, hologramIndex, backpropagate)
            % propagate Propagates a hologram using the Angular Spectrum method  
            %  
            %   hologramPostPropagation = propagate(obj, hologramIndex, willReplace) propagates  
            %   the hologram at the specified index using the Angular Spectrum method. If  
            %   the propagation distance is negative, the propagation is inverted. The result  
            %   is returned as a complex-valued hologram. Optionally, the propagated hologram  
            %   can replace the original one if the `willReplace` flag is true.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - hologramIndex: (scalar double) Index of the hologram to propagate.  
            %       - willReplace: (logical) If true, replaces the original hologram with the  
            %         propagated one.  
            %  
            %   Output:  
            %       - hologramPostPropagation: (M×N double) The propagated hologram image.
            if(nargin <3) 
                backpropagate = 0; 
            end

            kernel = obj.PropagationKernels(:,:,hologramIndex);

            if backpropagate
                FT_Uout = kernel.*(fft2(conj(obj.Holograms(:,:, hologramIndex))));
                hologramPostPropagation = conj(ifft2(FT_Uout));
            else
                FT_Uout = kernel.*fft2(obj.Holograms(:,:, hologramIndex));
                hologramPostPropagation = ifft2(FT_Uout);
            end
            
            hologramPostPropagation = double(hologramPostPropagation);
        end

        function propagatedField = propagateFor(obj, inputField, propagationDistance, wavelength)
            % Custom propagation using angular spectrum method based on QCI_Model parameters.
            waveNumber = 2 * pi / wavelength;
        
            fx = obj.XInverseCoordinates;
            fy = obj.YInverseCoordinates;
        
            % Meshgrid reuse
            [FX, FY] = meshgrid(fx, fy);
            if(propagationDistance < 0)
                phaseKernel = exp(-1i * propagationDistance * waveNumber * ...
                    sqrt(obj.MediumRefractiveIndex^2 - wavelength^2 * (FX.^2 + FY.^2)));
            
                inputFieldFFT = phaseKernel .* fftshift(fft2(fftshift(conj(inputField))));
                propagatedField = conj(fftshift(ifft2(ifftshift(inputFieldFFT))));
            else
                phaseKernel = exp(1i * propagationDistance * waveNumber * ...
                    sqrt(obj.MediumRefractiveIndex^2 - wavelength^2 * (FX.^2 + FY.^2)));
            
                inputFieldFFT = phaseKernel.*fftshift(fft2(fftshift(inputField)));
                propagatedField = fftshift((ifft2(ifftshift(inputFieldFFT))));
            end
        end

        function [fixedReferenceImages, fixedHolograms, correctedHeights, tforms] = threePointShiftAndScalingCorrection(obj, referenceImages, points)
            % threePointShiftAndScalingCorrection Corrects holograms using a three-point shift  
            % and scaling method  
            %  
            %   [fixedReferenceImages, fixedHolograms, correctedHeights, tforms] =  
            %   threePointShiftAndScalingCorrection(obj, referenceImages, points) applies a  
            %   three-point shift and scaling correction to a set of holograms and reference  
            %   images. The function computes a transformation based on matching points and  
            %   applies it to the holograms. The corrected heights are also calculated based  
            %   on the transformations.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - referenceImages: (M×N×P double) Set of reference images for correction.  
            %       - points: Cell array containing sets of corresponding points for each image.  
            %  
            %   Outputs:  
            %       - fixedReferenceImages: (M×N×P double) The corrected reference images.  
            %       - fixedHolograms: (M×N×P double) The corrected holograms.  
            %       - correctedHeights: (1×P double) The corrected propagation distances (heights).  
            %       - tforms: Cell array of transformation objects applied to the holograms. 
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
            % automaticalShiftAndScalingCorrection Automatically corrects holograms using  
            % feature-based shift and scaling alignment  
            %  
            %   [fixedReferenceImages, fixedHolograms, correctedHeights, tforms] =  
            %   automaticalShiftAndScalingCorrection(obj, referenceImages) automatically aligns  
            %   a series of reference images and holograms by detecting and matching features.  
            %   The transformations are estimated using a similarity transform, and the holograms  
            %   are corrected. The corrected heights are calculated based on the transformations.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - referenceImages: (M×N×P double) Set of reference images for alignment.  
            %  
            %   Outputs:  
            %       - fixedReferenceImages: (M×N×P double) The corrected reference images.  
            %       - fixedHolograms: (M×N×P double) The corrected holograms.  
            %       - correctedHeights: (1×P double) The corrected propagation distances (heights).  
            %       - tforms: Cell array of transformation objects applied to the holograms.  
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
            referenceImageFixed = preprocessHologram(referenceImages(:, :, end));
            fixedReferenceImages(:, :, end) = referenceImageFixed;
            fixedHolograms(:, :, end) = holograms(:, :, end);
            correctedHeights(end) = obj.PropagationDistances(end);
            
            % Detect features in the fixed reference image
            estimator = getAffineEstimator(referenceImageFixed, "KAZE");
            %estimator.params.MetricThreshold= 10000;
        
            % --- Loop through all images except the last one ---
            for k = 1:imageNumber-1
                % Current distorted reference image and hologram
                distortedImage = preprocessHologram(referenceImages(:, :, k));
                distortedHologram = holograms(:, :, k);
        
                tform = estimator.getAffineMatrix(distortedImage);
        
                % Step 4: Apply transformation to images and holograms
                Routput = imref2d(size(referenceImageFixed));
                fixedReferenceImages(:, :, k) = imwarp(distortedImage, tform, 'OutputView', Routput);
                fixedHolograms(:, :, k) = imwarp(distortedHologram, tform, 'OutputView', Routput);
                
                %showResult(referenceImageFixed, distortedImage, tform); %TEMP

                tforms{k} = tform;
                correctedHeights(k) = round(obj.PropagationDistances(k).*tform.T(1,1).^2);
        
                % Optional Debugging (uncomment for visualization)
                % figure; showMatchedFeatures(referenceImageFixed, distortedImage, ...
                %                            matchedFixed(inlierIdx), matchedDistorted(inlierIdx));
                % title(sprintf('Image %d Matching Points (Inliers Only)', k));
            end
        end

        function [fixedHolograms, rowShifts, columnShifts] = automaticalShiftCorrection(obj, holograms)
            % automaticalShiftCorrection Automatically corrects shifts in holograms  
            %  
            %   [fixedHolograms, rowShifts, columnShifts] = automaticalShiftCorrection(obj, holograms)  
            %   applies a phase correlation method to correct horizontal and vertical shifts  
            %   between holograms. The holograms are aligned by compensating for shifts in  
            %   rows and columns. The corrected holograms are returned along with the calculated  
            %   shifts for each hologram.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - holograms: (M×N×P double) The set of holograms to be aligned.  
            %  
            %   Outputs:  
            %       - fixedHolograms: (M×N×P double) The corrected holograms after alignment.  
            %       - rowShifts: (1×P double) The vertical shifts for each hologram.  
            %       - columnShifts: (1×P double) The horizontal shifts for each hologram.  
            if(isempty(holograms))
                error("No holograms to work with.");
            end

            hologramsFFT = fft2(holograms);
            hologramCount = size(holograms,3);

            usfactor = 10;
            rowShifts = zeros(1,hologramCount);
            columnShifts = zeros(1,hologramCount);

            for hologramIndex = 1:(hologramCount-1)
                [output, Greg] = obj.dftregistration(hologramsFFT(:,:,end), hologramsFFT(:,:,hologramIndex), usfactor);

                rowShifts(hologramIndex) = output(3);
                columnShifts(hologramIndex) = output(4);
                fixedHolograms(:,:,hologramIndex) = ifft2(Greg);
            end

            rowShifts(end) = 0;
            columnShifts(end) = 0;
            fixedHolograms(:,:,end+1) = holograms(:,:,end);
        end

        function shiftedHologram = shiftHologram(obj, hologram, shift_x, shift_y)
            % shiftHologram Shifts a hologram by specified amounts in the x and y directions
            %  
            %   shiftedHologram = shiftHologram(obj, hologram, shift_x, shift_y)  
            %   applies a phase shift in the frequency domain to shift a given hologram by  
            %   the specified amounts in the horizontal (x) and vertical (y) directions.  
            %  
            %   Inputs:  
            %       - obj: Instance of QCI_Model.  
            %       - hologram: (M×N double) The hologram to be shifted.  
            %       - shift_x: (double) The shift in the x-direction (horizontal).  
            %       - shift_y: (double) The shift in the y-direction (vertical).  
            %  
            %   Outputs:  
            %       - shiftedHologram: (M×N double) The shifted hologram. 

            [Ny, Nx] = size(hologram); % Hologram size
             
            % Siatki przestrzeni częstotliwości
            fx = [0:Nx/2-1, -Nx/2:-1] / Nx; % Normalize x frequencies
            fy = [0:Ny/2-1, -Ny/2:-1] / Ny; % Normalize y frequencies
            [Fx, Fy] = meshgrid(fx, fy);
             
            % Phase shift in frequency domain
            phase_shift = exp(-1i * 2 * pi * (Fx * shift_x + Fy * shift_y));
             
            % Image shift
            shiftedHologram = ifft2(fft2(hologram) .* phase_shift);
        end

        function reconstruction = IGA(obj, iter, sigma)
            % Iterative Gabor Averagin (IGA) - method for phase retrieval from multiple
            % in-line holograms collected with different defocus and/or wavelength.
            % Algorithm was designed to work with low signal-to-noise-ratio data
            %
            % Inputs:
            %   iter - number of iterations. Default - iter = 5;
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
        
        function reconstruction = GerchbergSaxton(obj, iterationCount, gaussianSigma)
            % Optimized Gerchberg-Saxton reconstruction using local variables and custom propagation.
        
            hologramsOriginal = obj.Holograms;
        
            amplitudeImages = sqrt(hologramsOriginal);
            numberOfHolograms = size(hologramsOriginal, 3);
        
            inputField = amplitudeImages;
        
            for iter = 1:iterationCount
                for idx = 1:(numberOfHolograms - 1)
                    reconstruction = propagateOptimized(inputField(:,:, idx), obj.PropagationKernels(:,:,idx));

                    fieldAmplitude = abs(reconstruction);
                    if gaussianSigma > 0
                        fieldAmplitudeGauss = imgaussfilt(fieldAmplitude, gaussianSigma);
                        fieldAmplitude(fieldAmplitude>fieldAmplitudeGauss) = fieldAmplitudeGauss(fieldAmplitude>fieldAmplitudeGauss);
                    end

                    phase = angle(reconstruction);
                    
                    reconstruction = fieldAmplitude.*exp(1i*phase);

                    inputField(:, :, idx+1) = propagateOptimized(reconstruction, obj.PropagationKernels(:,:,idx+1), true);
                    inputField(:, :, idx+1) = amplitudeImages(:,:, idx+1).*inputField(:,:, idx+1)./abs(inputField(:, :, idx+1));
                end
                
                reconstruction = propagateOptimized(inputField(:,:,end), obj.PropagationKernels(:,:,end));

                inputField(:, :, 1) = propagateOptimized(reconstruction, obj.PropagationKernels(:,:,1), true);

                inputField(:, :, 1) = amplitudeImages(:, :, 1).*inputField(:, :, 1)./abs(inputField(:,:,1));
            end
                    
            % Final propagation to object plane
            reconstruction = propagateOptimized(inputField(:,:, 1), obj.PropagationKernels(:,:, 1));
            reconstruction = reconstruction.^2;
        end

        function reconstruction = GaborAveraging(obj)
            % GaborAveraging - Gabor averaging in-line holography reconstruction
            %
            %   reconstruction = GaborAveraging(obj)
            %   This method performs Gabor averaging, a technique used to reconstruct 
            %   an object field from multiple holograms by backpropagating each hologram 
            %   to the object plane and averaging the results.
            %
            %   The result is a reconstruction that incorporates the phase and amplitude 
            %   of the optical field from all holograms.
            %
            %   Outputs:
            %       - reconstruction: (M×N double) The final averaged object field reconstruction
            %                         after applying the Gabor averaging technique.
            reconstruction = zeros(size(obj.Holograms(:,:,1)));
            hologramCount = size(obj.Holograms, 3);

            for nn = 1:hologramCount
                % Backpropagate each hologram to object plane
                intermediateReconstruction = obj.propagate(nn, true);
                % Add the propagation result to the R_GA (and rescale phase to
                % match the 1st wavelength)
                reconstruction = reconstruction + sqrt(abs(intermediateReconstruction)) .* exp(1i.*angle(intermediateReconstruction) * obj.Wavelengths(nn)/obj.Wavelengths(1));
            end

            % Divide by the number of holograms
            reconstruction = reconstruction/hologramCount;

            obj.PropagationDistances = -obj.PropagationDistances;
        end

        % DarkFocus version 3
        function [bestFocusZ, focusCurves, DarkVolume] = DarkFocus(obj, hologramIndex, range, inspectedROI, darkVolumeROI, isAmplitude)
            % DarkFocus - Optimizing focus by evaluating the variance of gradient magnitudes in the dark volume
            %
            %   [bestFocusZ, focusCurves, DarkVolume] = DarkFocus(obj, hologramIndex, range, inspectedROI, darkVolumeROI)
            %
            %   This method evaluates the sharpness of a hologram over a given range of focus 
            %   by calculating the variance of the gradient magnitudes in a selected dark volume region.
            %   The goal is to identify the best focus by finding the range that maximizes focus sharpness.
            %
            %   Inputs:
            %       - hologramIndex: (int) The index of the hologram to process.
            %       - range: (1xN array) The range of propagation distances (Z-steps) to inspect for focus.
            %       - inspectedROI: (1x4 array) The region of interest (ROI) within the hologram for inspection.
            %       - darkVolumeROI: (1x4 array) The ROI for calculating the dark volume (default is the whole inspectedROI).
            %
            %   Outputs:
            %       - bestFocusZ: (double) The propagation distance (Z-coordinate) corresponding to the best focus.
            %       - focusCurves: (1xN array) The calculated focus curves showing sharpness over the range.
            %       - DarkVolume: (MxN array) The amplitude of the dark volume at the best focus.
            %
            %   The function uses a model to propagate the hologram through the range and computes
            %   the variance of gradient magnitudes at each step to find the best focus.

            background = imgaussfilt(obj.Holograms(:,:,hologramIndex), 30);
            darkHologram = obj.Holograms(:,:,hologramIndex) - background;

            YIndexes = (inspectedROI(2):inspectedROI(2)+inspectedROI(4) - 1) + 1;
            XIndexes = (inspectedROI(1):inspectedROI(1)+inspectedROI(3) - 1) + 1;

            if(darkVolumeROI == 0)
                darkVolumeROI = [1 1 inspectedROI(3) inspectedROI(4)];
            else
                darkVolumeROI(2) = darkVolumeROI(2) - inspectedROI(2) + 1;
                darkVolumeROI(1) = darkVolumeROI(1) - inspectedROI(1) + 1;
            end

            darkVolumeYIndexes = (darkVolumeROI(2):(darkVolumeROI(2)+darkVolumeROI(4)) - 1);
            darkVolumeXIndexes = (darkVolumeROI(1):(darkVolumeROI(1)+darkVolumeROI(3)) - 1);
            
            darkHologram = darkHologram(YIndexes,XIndexes);

            darkModel = QCI_Model("darkHologram",darkHologram,obj.Wavelengths(hologramIndex),obj.PropagationDistances(hologramIndex),obj.CameraPixelSize);
            DarkFocus = zeros(1,length(range));
            DarkVolume = zeros(length(darkVolumeYIndexes), length(darkVolumeXIndexes));

            for rangeIndex = 1:length(range)
                darkModel.setPropagationDistances(range(rangeIndex),1);
                Obj1 = darkModel.propagate(1);
                % DarkFocus
                if isAmplitude == 1
                    img = abs(Obj1);
                else
                    img = angle(Obj1);
                end
                DarkVolume(:,:,rangeIndex) = img(darkVolumeYIndexes, darkVolumeXIndexes); %angle(Obj1);%
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

        function [output, Greg] = dftregistration(obj, buf1ft,buf2ft,usfac)
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
                CC = ifft2(obj.FTpad(buf1ft.*conj(buf2ft),[2*nr,2*nc]));
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
                    CC = conj(obj.dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
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

        function out=dftups(obj, in,nor,noc,usfac,roff,coff)
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
        
        function [ imFTout ] = FTpad(obj, imFT,outsize)
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

function propagatedField = propagateOptimized(inputField, kernel, isBackpropagated)
    if nargin < 3
        isBackpropagated = 0;
    end
    
    if isBackpropagated
        fieldFFT = kernel.*(fft2(conj(inputField)));
        propagatedField = conj(ifft2(fieldFFT));
    else
        fieldFFT = kernel.*fft2(inputField);
        propagatedField = ifft2(fieldFFT);
    end
end

function imageProcessed = preprocessHologram(image)
    imageProcessed = (angle(image) + pi)/(2*pi);
    imageProcessed = imgradient(imageProcessed, 'Sobel');
end


function showResult(imageOrigin, imageDistorted, tform) %TEMP
    outputView = imref2d(size(imageOrigin));
    warpedImage = imwarp(imageDistorted, tform, 'OutputView', outputView);

    figure()
    imshowpair(imageOrigin, warpedImage, 'falsecolor');
    colormap('jet');
    title('Result of affine transformation');
end