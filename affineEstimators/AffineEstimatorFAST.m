classdef AffineEstimatorFAST < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorFAST(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.MinQuality = 0.1;
            obj.params.MinContrast = 0.2;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectFASTFeatures(image, ...
                                        'MinQuality', obj.params.MinQuality, ...
                                        'MinContrast', obj.params.MinContrast);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
