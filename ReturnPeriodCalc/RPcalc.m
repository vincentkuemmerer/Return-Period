function [R,HRP,Lower_RP,Upper_RP,icStorm,info,TT_POT,ST,ID,IT] = RPcalc(Data,WinterMonths,TempRes,MSD,THperc,CF_alpha)
    warning('off','all')    
    
    % Data: timetable(Time,Hs)
    % WinterMonths = [12,1,2,3]; % Winter season DJFM, Northern Hemisphere or [7,8,9,10]; % Winter season JASO, Southern Hemisphere
    % TempRes: Temporal Resulotion (Sampling interval)
    % MSD: Minimum Storm Duration (6hours)
    % THperc: Percentage of storm threshold (ST)
    % CF_alpha: Confidence interval for return period calculation

    % --------------------------------------------------- Get Winter Season    

    DataWinter = Data(ismember(month(Data.Time),WinterMonths),:);    
    HsW = DataWinter.Hs; % Significant Wave Height in Winter

    % --------------------------------------------------- Storm Identification 

    ST = prctile(Data.Hs,THperc); % Threshold defition by the 95th percentile of the entire data
    IT = mean(DataWinter.Hs,'omitnan');
    [~,~,Z] = ExtremalIndex(HsW,ST); % Extremal index inference    
    ID = RunLength(HsW,ST,1:750,Z,TempRes); % Run length (independence criterion for consecutive storms) inference

    if ~isnan(Z) && ~isempty(HsW)
    [icStorm,~,info] = StormIdentification(HsW,TempRes,ST,ID,IT,MSD); % Storm Identification

    if ~isempty(icStorm) && numel(icStorm) > 1

        % --------------------------------------------------- Peak Over Threshold (POT) method

        POT = zeros(1,length(icStorm));
        l = zeros(1,length(icStorm));
        timePOT = NaT(length(icStorm),1);
        for i = 1:length(icStorm)
            [POT(i),l(i)] = max(HsW(icStorm{i})); % Get the max value of the individual storms
            timePOT(i) = Data.Time(icStorm{i}(l(i))); % Get the time for the max POT value of each storm
        end
        POT = POT'-ST; % Substract the TH level from the extreme value data set
        TT_POT = timetable(timePOT,POT); 
        
        % Maximum Likelihood Estimation of GPD parameters
        % Significance level for the confidence interval pci of parameter estimates, specified as a scalar in the range (0,1). 
        % The confidence level of pci is 100(1–Alpha)%. The default is 0.05 for 95% confidence. 
        % Example: 'Alpha',0.01 specifies the confidence level as 99%.

        [MLE,MLE_CFI] = mle(POT,'distribution','Generalized Pareto','theta',0,'Alpha',CF_alpha); % 50% confidence level
        shape = MLE(1); 
        scale = MLE(2);

        % Confidence intervals
        lower_shape = MLE_CFI(1,1); 
        lower_scale = MLE_CFI(1,2);
        upper_shape = MLE_CFI(2,1);
        upper_scale = MLE_CFI(2,2);
        
        % Return Period - Buoy Data
        Z2 = length(POT);
        lambda = range(year(Data.Time));
        R = [1:10,20:10:100];
        for i = 1:length(R)
            HRP(i) = ST + (scale/shape)*((((R(i)*Z2)/lambda)^shape)-1);
            Lower_RP(i) =  ST + (lower_scale/lower_shape)*((((R(i)*Z2)/lambda)^lower_shape)-1);
            Upper_RP(i) =  ST + (upper_scale/upper_shape)*((((R(i)*Z2)/lambda)^upper_shape)-1);
        end

    else
        R=0;HRP=0;Lower_RP=0;Upper_RP=0;TT_POT=0;
    end
    else
        icStorm=0;info=0;R=0;HRP=0;Lower_RP=0;Upper_RP=0;TT_POT=0;
    end
end