classdef AffineEstimatorSURF < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorSURF(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.MetricThreshold = 1000;
            obj.params.NumOctaves = 3;
            obj.params.NumScaleLevels = 4;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectSURFFeatures(image, ...
                                        'MetricThreshold', obj.params.MetricThreshold, ...
                                        'NumOctaves', obj.params.NumOctaves, ...
                                        'NumScaleLevels', obj.params.NumScaleLevels);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
