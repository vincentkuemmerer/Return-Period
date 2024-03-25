function [r,z] = RunLength(Hs,Threshold,r_estimate,Z,Resolution)
    W = Hs>Threshold;
    Res = Resolution;
    %r_estimate = 1:750; % Default
    z = zeros(1,length(r_estimate));
    for k = 1:length(r_estimate)
        r = r_estimate(k);
        % disp(['r=',num2str(ri)])
        for i = 2:length(W)-r
            if (W(i-1) == 1) && (W(i) == 0)
                if sum(W(i:i+r-1)) == 0
                    z(k) = z(k)+1;
                end
            end
        end
         if z(k) <= round(Z)
            z = z(k);
            r = r_estimate(k)*Res;
            % disp(['z=',num2str(z(k))]);
            break
        end
    end
end