function [Tmid, Nend] = EggHalftime(cumEgg)
NumWorm = size(cumEgg, 1);
NumT = size(cumEgg, 2);
Tmid = nan(NumWorm, 1);
Nend = nan(NumWorm, 1);
T = (0: (NumT-1))/30;
for i = 1 : NumWorm
    currcumEgg = cumEgg(i, :);
    currmidEgg = repmat(currcumEgg(end)/2, [1 NumT]);
    vcumEgg = [T; currcumEgg];
    vmidEgg = [T; currmidEgg];
    currP = InterX(vcumEgg, vmidEgg);
    Nend(i) = currcumEgg(end)/T(end);
    if ~isempty(currP)
        Tmid(i) = currP(1);
    end
end
% Tmid(isnan(Tmid)) = [];
end

