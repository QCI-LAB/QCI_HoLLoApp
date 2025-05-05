classdef AffineEstimatorORB < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorORB(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.ScaleFactor = 1.2;
            obj.params.NumLevels = 8;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectORBFeatures(image, ...
                                        'ScaleFactor', obj.params.ScaleFactor, ...
                                        'NumLevels', obj.params.NumLevels);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
