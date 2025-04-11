classdef AffineEstimatorKAZE < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorKAZE(image_base)
            obj = obj@AffineEstimatorBase([]);

            obj.rescale_factor_base = 0.5;
            obj.rescale_factor_new = 0.5;

            obj.params.Threshold = 0.0001;
            obj.params.NumOctaves = 3;
            obj.params.NumScaleLevels = 4;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectKAZEFeatures(image, ...
                                        'Threshold', obj.params.Threshold, ...
                                        'NumOctaves', obj.params.NumOctaves, ...
                                        'NumScaleLevels', obj.params.NumScaleLevels);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
