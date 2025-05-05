function detector = getAffineEstimator(image_base, name_method)
    % getAffineEstimator Estimates an affine transformation detector based on the selected method.
    %
    % Arguments:
    %   image_base (required) : The image to be used for affine estimation.
    %   name_method (required) : The method used for affine transformation estimation. 
    %                            Valid options are 'SURF', 'BRISK', 'FAST', 'ORB', 
    %                            'MinEigen', 'Harris', 'KAZE', and 'SIFT'.
    %
    % Returns:
    %   detector : A detector object created using the specified method, or an empty array 
    %              if an invalid method is selected. If an invalid method is provided, 
    %              a warning is issued.
    %
    % Example:
    %   detector = getAffineEstimator(image, 'SURF');
    
    arguments
        image_base % Automatically detects the type (e.g., image)
        name_method {mustBeMember(name_method, {'SURF', 'BRISK', 'FAST', 'ORB', 'MinEigen', 'Harris', 'KAZE', 'SIFT'})} % Argument validation
    end    

    switch name_method
        case 'SURF'
            detector = AffineEstimatorSURF(image_base);
        case 'BRISK'
            detector = AffineEstimatorBRISK(image_base);
        case 'FAST'
            detector = AffineEstimatorFAST(image_base);
        case 'ORB'
            detector = AffineEstimatorORB(image_base);
        case 'MinEigen'
            detector = AffineEstimatorMinEigen(image_base);
        case 'Harris'
            detector = AffineEstimatorHarris(image_base);
        case 'KAZE'
            detector = AffineEstimatorKAZE(image_base);
        case 'SIFT'
            detector = AffineEstimatorSIFT(image_base);
        otherwise
            warning('Unknown method "%s", returning empty detector.', name_method);
            detector = [];
    end
end
