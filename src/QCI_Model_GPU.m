classdef QCI_Model_GPU < QCI_Model
    %QCI_MODEL_GPU Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods
        function obj = QCI_Model_GPU(names, holograms, wavelengths, propagationDistances, pixelSize)
           obj = obj@QCI_Model(names, holograms, wavelengths, propagationDistances, pixelSize);
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

        function addHologram(obj, name, hologram, wavelength, propagationDistances)
            obj.Names(end + 1) = name;
            obj.Holograms(:, :, end + 1) = gpuArray(hologram);
            obj.Wavelengths(end + 1) = wavelength;
            obj.PropagationDistances(end + 1) = propagationDistances;
            obj.WaveNumbers(end + 1) = 2*pi./obj.Wavelengths(end);
            
            obj.HologramsFFT(:,:,end + 1) = fft2(obj.Holograms(:,:,end));
            [ySize, xSize] = size(obj.Holograms(:, :, 1));
            obj.FreeSpaceImpulseMatrixes(:,:,end + 1) = fftshift(obj.WaveNumbers(end) *...
                sqrt(obj.MediumRefractiveIndex^2 - obj.Wavelengths(end)^2 *...
                (ones(ySize,1)*(obj.XInverseCoordinates.^2) + (obj.YInverseCoordinates'.^2)*ones(1,xSize))));
        end

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

            kernel = gpuArray(obj.PropagationKernels(:,:,hologramIndex));

            if backpropagate
                FT_Uout = kernel.*(fft2(conj(obj.Holograms(:,:, hologramIndex))));
                hologramPostPropagation = conj(ifft2(FT_Uout));
            else
                FT_Uout = kernel.*fft2(obj.Holograms(:,:, hologramIndex));
                hologramPostPropagation = ifft2(FT_Uout);
            end
            
            hologramPostPropagation = double(gather(hologramPostPropagation));
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
            if(numberOfHolograms > 1)
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
            else
                for iter = 1:iterationCount
                    reconstruction = propagateOptimized(inputField(:,:, 1), obj.PropagationKernels(:,:,1));

                    fieldAmplitude = abs(reconstruction);
                    phase = angle(reconstruction);
                    
                    if gaussianSigma > 0
                        fieldAmplitudeGauss = imgaussfilt(fieldAmplitude, gaussianSigma);
                        fieldAmplitude(fieldAmplitude>fieldAmplitudeGauss) = fieldAmplitudeGauss(fieldAmplitude>fieldAmplitudeGauss);
                    end

                    reconstruction = fieldAmplitude.*exp(1i*phase);

                    inputField(:, :, 1) = propagateOptimized(reconstruction, obj.PropagationKernels(:,:,1), true);
                    inputField(:, :, 1) = amplitudeImages(:,:, 1).*inputField(:,:, 1)./abs(inputField(:, :, 1));
                end
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