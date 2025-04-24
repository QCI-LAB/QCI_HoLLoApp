classdef AffineEstimatorSIFT < AffineEstimatorBase
    properties
        params
    end

    methods
        function obj = AffineEstimatorSIFT(image_base)
            obj = obj@AffineEstimatorBase([]);  
            
            obj.params.ContrastThreshold = 0.0133;
            obj.params.EdgeThreshold = 10.0;
            obj.params.NumLayersInOctave = 3;
            obj.params.Sigma = 1.6;

            if ~isempty(image_base)
                obj.setNewImage(image_base); % teraz już można bezpiecznie
            end
        end

        function [features_new, validPoints_new] = detectAndExtractFeatures(obj, image)
            points = detectSIFTFeatures(image, ...
                                        'ContrastThreshold', obj.params.ContrastThreshold, ...
                                        'EdgeThreshold', obj.params.EdgeThreshold, ...
                                        'NumLayersInOctave', obj.params.NumLayersInOctave, ...
                                        'Sigma', obj.params.Sigma);
            [features_new, validPoints_new] = extractFeatures(image, points);
        end
    end
end
