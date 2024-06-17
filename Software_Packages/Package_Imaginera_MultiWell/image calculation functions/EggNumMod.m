function EggNum_mod = EggNumMod(EggNum)
EggNum_mod = zeros(size(EggNum));
NumW = size(EggNum, 1);
NumT = size(EggNum, 2);
for i = 1 : NumW
    V = EggNum(i, :);
    for j = 2 : NumT
        if V(j)<V(j-1)
            V(j) = V(j-1);
        end
    end
    EggNum_mod(i, :) = V;
end
end

