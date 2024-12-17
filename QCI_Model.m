classdef QCI_Model < handle
    %QCI_MODEL is a bundle used for hologram creation and manipulation
    %   Detailed explanation goes here
    
    properties (Access = private)
        %Initial Data
        Names                   (1,:) string = []
        ImagesRaw               (:,:,:) double                              % Hologram image (matrix)
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
        ImagesFFT               (:,:,:)                          = []
        ImagesPostPropagation   (:,:,:)                          = []
        
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
                obj.ImagesRaw(:,:,i) = double(cell2mat(images{i}));
            end
            obj.Wavelengths = wavelengths;
            obj.ZCoordinates = distances;
            obj.CameraPixelSize = pixelSize;
            [ySize, xSize] = size(obj.ImagesRaw(:,:,1));
            
            % Deriving data
            obj.YCoordinates = (0:(ySize-1))*obj.CameraPixelSize;
            obj.XCoordinates = (0:(xSize-1))*obj.CameraPixelSize;
            inverseStepY = 1/ySize/pixelSize;
            inverseStepX = 1/xSize/pixelSize;
            YInverseCoordinates =  -ySize/2*inverseStepY:inverseStepY:(ySize/2 - 1)*inverseStepY;
            XInverseCoordinates =  -xSize/2*inverseStepX:inverseStepX:(xSize/2 - 1)*inverseStepX;
            obj.WaveNumbers = 2*pi./obj.Wavelengths;
            for i = 1:length(wavelengths)
                obj.ImagesFFT(:,:,i) = fft2(obj.ImagesRaw(:,:,i));
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
            image = uint8(obj.ImagesRaw(:,:, hologramIndex));
            return
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

        function imagePostPropagation = propagate(obj,z,imageIndex)
            % Propagates the image with Angular Scaling method
            kernelExponent = obj.FreeSpaceImpulseMatrixes(:,:,imageIndex)*z;
            kernelExponent = kernelExponent - kernelExponent(1,1);
            
            if(z<0)
                kernel = exp(-1i*kernelExponent);
                propagetedFFT = kernel.*obj.ImagesFFT(:,:, imageIndex);
                obj.ImagesPostPropagation(:,:, imageIndex) = conj(ifft2(propagetedFFT));
                imagePostPropagation = obj.ImagesPostPropagation(:,:, imageIndex);
            else
                kernel = exp(1i*kernelExponent);
                propagetedFFT = kernel.*obj.ImagesFFT(:,:, imageIndex);
                obj.ImagesPostPropagation(:,:, imageIndex) = ifft2(propagetedFFT);
                imagePostPropagation = obj.ImagesPostPropagation(:,:, imageIndex);
            end
        end
    end
end

