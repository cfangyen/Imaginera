function D = Str2Mat(A)
%STR2MAT convert string to matrix
D = reshape(str2double(regexp(A,'\d*','match')),2,[])';

end

