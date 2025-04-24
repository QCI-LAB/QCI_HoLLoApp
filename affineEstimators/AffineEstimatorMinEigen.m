classdef AffineEstimatorMinEigen < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorMinEigen(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.MinQuality = 0.1;
            obj.params.FilterSize = 5;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectMinEigenFeatures(image, ...
                                        'MinQuality', obj.params.MinQuality, ...
                                        'FilterSize', obj.params.FilterSize);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
