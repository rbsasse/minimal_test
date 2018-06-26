function [hash] = util_md5_external(string)
% generate a md5 hash from given text string using external command
% note: windows probably won't have it
% matlab doesn't seem to have any lightweight solution

archstr = computer('arch');

switch archstr
    
    case 'win64'
        
        warning('Windows, no md5 command.')
        command = sprintf('echo -n %s', string);
        
    case 'maci64'
        
        command = sprintf('echo -n %s | md5', string);
        
    case 'glnxa64'
        
        command = sprintf('echo -n %s | md5sum', string);
        
    otherwise
        
        warning('Unexpected OS type.')
        command = sprintf('echo -n %s', string);
        
end

[status,cmdout] = system(command);
if status ~= 0
    error('md5 command failed')
end

hash = regexprep(cmdout,'\r\n|\n|\r| -| ',''); % remove trailing chars
end
