classdef AffineEstimatorHarris < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorHarris(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.MinQuality = 0.1;
            obj.params.FilterSize = 5;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectHarrisFeatures(image, ...
                                        'MinQuality', obj.params.MinQuality, ...
                                        'FilterSize', obj.params.FilterSize);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
