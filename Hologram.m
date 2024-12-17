classdef Hologram < handle
    properties
        Image (:,:)                                                         % Hologram image (matrix)
        Wavelength (1,1) double {mustBeNonnegative} = 0.550                 % Wavelength of light used [micro meters]
        Distance   (1,1) double {mustBeNonnegative} = 0                     % Distance from CCD to sample [micro meters]
        MediumRefractiveIndex (1,1) double = 1                              % Refractive index of the medium in which the hologram was taken (assume air)
        CameraPixelSize  (1,1) double = 1                                      % Size of a pixel on the CCD camera used to obtain the hologram
    end
    
    methods
        %% Constructor
        function obj = Hologram(image, wavelength, distance, pixelSize)
            if nargin == 0
                obj.Image = [];
                obj.Wavelength = 0;
                obj.Distance = 0;
                obj.CameraPixelSize = 0;
                return
            end
            obj.Image = image;
            obj.Wavelength = wavelength;
            obj.Distance = distance;
            obj.CameraPixelSize = pixelSize;
        end
        
        %% Focus Adjustment
        function obj = focus(obj)
            disp("Adjusting focus of the hologram...");
            obj.Image = focusHologram(obj.Image, obj.Distance, obj.Wavelength);
        end
        
        %% Remove Shift
        function obj = removeShift(obj)
            disp("Removing shift from the hologram...");
            obj.Image = removeHologramShift(obj.Image);
        end
        
        %% Scale Hologram
        function obj = scale(obj, referenceDistance)
            disp("Scaling the hologram to a reference distance...");
            obj.Image = scaleHologram(obj.Image, obj.Distance, referenceDistance);
        end
    end
end

