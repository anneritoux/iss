%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss_OMP; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values
classdef iss_OMP_ConstantBackground_WeightDotProduct < iss_OMP
    
    properties 
        %To stop blow up when using weights.
        ompWeightShift = 0.01;
        
        %1.0 means each round treated equally, The lower the value,
        %the greater the influence of the round with the largest residual.
        ompWeightPower = 0.9;     
                
        %Need to normalise bled codes so rounds where gene efficiencies low
        %don't contribute when working out which gene to add next. 
        %Rounds With GeneEfficiencies less than ompNormBledCodeShift won't
        %contribute very much. 
        %ompNormBledCodeScale determines the sharpness of the drop off in
        %influence i.e. larger value means low gene efficiency rounds
        %contribute less. 
        ompNormBledCodeShift = 0.5;
        ompNormBledCodeScale = 7;
        ompNormBledCodeUnbledBoost = 1;
        
    end
end
