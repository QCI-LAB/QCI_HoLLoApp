classdef AffineEstimatorBRISK < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorBRISK(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.MinContrast = 0.2;
            obj.params.MinQuality = 0.1;
            obj.params.NumOctaves = 4;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectBRISKFeatures(image, ...
                                        'MinContrast', obj.params.MinContrast, ...
                                        'MinQUality', obj.params.MinQuality, ...
                                        'NumOctaves', obj.params.NumOctaves);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
