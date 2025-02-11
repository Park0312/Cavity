function [sortedStrings, index] = natsortfiles(strings)
    [~, index] = sort(str2double(regexp(strings, '\d+\.?\d*', 'match', 'once')));
    sortedStrings = strings(index);
end