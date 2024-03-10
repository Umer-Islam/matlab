%function file for A1Q4
function [unsortedVector, sortedVector] = func_A1Q4(m1)
    unsortedVector = m1(:);

    sortedVector = sort(unsortedVector, 'descend');
end
