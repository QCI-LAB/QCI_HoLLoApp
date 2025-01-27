classdef QCI_Model_GPU < QCI_Model
    %QCI_MODEL_GPU Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods
        function obj = QCI_Model_GPU(names, holograms, wavelengths, zCoordinates, pixelSize)
           holograms = gpuArray(holograms);
           obj@QCI_Model(names, holograms, wavelengths, zCoordinates, pixelSize);

           obj.HologramsFFT = gpuArray(obj.HologramsFFT);
           obj.FreeSpaceImpulseMatrixes = gpuArray(obj.FreeSpaceImpulseMatrixes);
        end

        function holograms = getHolograms(obj, indexes)
            if(indexes == 0)
                holograms = gather(obj.Holograms);
            else
                holograms = gather(obj.Holograms(:,:, indexes));
            end
        end

        function holograms = setHolograms(obj, holograms, indexes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            if(size(holograms, 3) == size(obj.Holograms,3))
                    obj.Holograms = gpuArray(double(holograms));
                    obj.HologramsFFT = fft2(holograms);
            elseif(size(holograms, 3) == size(indexes))
                obj.Holograms(:,:,indexes) = gpuArray(double(holograms));
                obj.HologramsFFT(:,:,indexes) = fft2(holograms);
            else
                error("Dimension of indexes array does not match the amount of given holograms.");
            end
        end

        function addHologram(obj, name, hologram, wavelength, zCoordinate)
            obj.Names(end + 1) = name;
            obj.Holograms(:, :, end + 1) = gpuArray(hologram);
            obj.Wavelengths(end + 1) = wavelength;
            obj.ZCoordinates(end + 1) = zCoordinate;
            obj.WaveNumbers(end + 1) = 2*pi./obj.Wavelengths(end);
            
            obj.HologramsFFT(:,:,end + 1) = fft2(obj.Holograms(:,:,end));
            [ySize, xSize] = size(obj.Holograms(:, :, 1));
            obj.FreeSpaceImpulseMatrixes(:,:,end + 1) = fftshift(obj.WaveNumbers(end) *...
                sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(end)^2 *...
                (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
        end
    end
end