classdef AffineEstimatorBase < handle
    % AffineEstimatorBase
    % 
    % This class provides an abstract base for estimating affine transformations
    % between a reference image and a new image. It performs feature detection,
    % extraction, matching, and transformation estimation. The final transformation 
    % can be rescaled according to different rescaling factors.
    properties (Access = protected)
        image_base = []
        image_base_resized = []
        features_base = []
        validPoints_base = []
    end

    properties 
        rescale_factor_base = 1
        rescale_factor_new = 1
        show_result = false
        fig1_num = 367
        fig2_num = 368
    end

    methods (Abstract)
        [features_new, validPoints_new] = detectAndExtractFeatures(obj, image);
    end
    
    methods
        function obj = AffineEstimatorBase(image_base)
            % Constructor method that initializes the base image if provided.
            if ~isempty(image_base)
                obj.setNewImage(image_base);
            end
        end

        function tform_final = getAffineMatrix(obj, image)
            % Computes the affine transformation matrix between the base image and the
            % new image. Optionally displays the results.
            image_resized = imresize(image, obj.rescale_factor_new);
            [features_new, validPoints_new] = obj.detectAndExtractFeatures(image_resized);
        
            [tform, matchedPoints_base, matchedPoints_new] = getTform(obj, features_new, validPoints_new);
        
            tform_final = obj.rescaleTform(tform);
            
            if obj.show_result
                obj.showMatchedPoints(image_resized, matchedPoints_base, matchedPoints_new);
                obj.showBlendedImages(image, tform_final);
            end
        end

        
        function setNewImage(obj, image_base)
            % Updates the base image, resizes it, and extracts its features.
            obj.image_base = image_base;
            obj.image_base_resized = imresize(image_base, obj.rescale_factor_base);

            [obj.features_base, obj.validPoints_base] = obj.detectAndExtractFeatures(obj.image_base_resized);
        end

        function set.rescale_factor_base(obj, value)
            % Setter for rescale_factor_base. Resizes the base image if necessary.
            obj.rescale_factor_base = value;
            if ~isempty(obj.image_base)
                obj.setNewImage(obj.image_base); 
            end
        end
    end
    
    methods (Access = protected)
        function [tform, matchedPoints_base, matchedPoints_new] = getTform(obj, features_new, validPoints_new)
            % Estimates the geometric transformation between the base image and the new image.
            indexPairs = matchFeatures(obj.features_base, features_new);
        
            matchedPoints_base = obj.validPoints_base(indexPairs(:, 1));
            matchedPoints_new = validPoints_new(indexPairs(:, 2));
        
            [tform, ~, ~] = estimateGeometricTransform2D(... 
                matchedPoints_new, matchedPoints_base, 'affine');
        end

        function showMatchedPoints(obj, image_resized, matchedPoints_base, matchedPoints_new)
            % Displays the matched points between the base and new images.
            figure(obj.fig1_num);
            showMatchedFeatures(image_resized, obj.image_base_resized, matchedPoints_new, matchedPoints_base, 'montage');
            title('Matched Points');
        end
        
        function showBlendedImages(obj, image, tform)
            % Displays the result of the affine transformation by blending the images.
            outputView = imref2d(size(obj.image_base));
            warpedImage = imwarp(image, tform, 'OutputView', outputView);
            
            figure(obj.fig2_num);
            imshowpair(obj.image_base, warpedImage, 'falsecolor');
            colormap('jet');
            title('Result of affine transformation');
        end

        function tform_final = rescaleTform(obj, tform)
            % Rescales the affine transformation matrix based on the new rescaling factor.
            scale_ratio = obj.rescale_factor_new / obj.rescale_factor_base;
            matrix = tform.T;

            % Rescale the linear part of the transformation matrix (rotation, scaling, shear)
            matrix_final = zeros(3, 3);
            matrix_final(3,3) = 1;
            matrix_final(1:2, 1:2) = matrix(1:2, 1:2) * scale_ratio;

            % Rescale the translation part of the transformation matrix
            matrix_final(3, 1:2) = matrix(3, 1:2) * (1 / obj.rescale_factor_base);

            % Create a new affine transformation object with the adjusted matrix
            tform_final = affine2d(matrix_final);
        end
    end
end
