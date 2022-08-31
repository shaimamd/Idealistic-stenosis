% Get all PDF files in the current folder
files = dir('*.f*');
% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
      % Convert to number
%       num = str2double(f);
%       if ~isnan(num)
          % If numeric, rename
          movefile(files(id).name, sprintf('stspipe0.f%05d', id));
%       end
end