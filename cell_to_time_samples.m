function A = cell_to_time_samples(xcells,idx)
% For a set of time indexes, create a matrix A that samples the rows of
% xcells assuming that the order of the cells corresponds to successive
% samples in time
% Inputs:
% - xcells = trials arranged into individual cells
% - idx = array of time points to include (default: all indexes)
% Outputs:
% - A = concatenated array or matrix of all included time points
% Nate Zuk (2018)

if size(xcells,1)==1, % if it's a row array
    xcells = xcells'; % turn it into a column array
end

allidx = cellfun(@(x) size(x,1),xcells); % length of time (indexes) for each cell
if nargin<2, idx = sum(allidx); end % if idx isn't specified, use all time (indexes)
ncols = size(xcells{1},2); % number of columns
A = zeros(length(idx),ncols); % preallocate matrix

% Bin idx by cell
sumallidx = [0; cumsum(allidx)];
for cc = 1:length(sumallidx)-1,
    idx_in_cell = idx>sumallidx(cc)&idx<=sumallidx(cc+1); % get indexes in the cell
    ntm = idx(idx_in_cell)-sumallidx(cc); % get the appropriate indexes within the cell
    A(idx_in_cell,:) = xcells{cc}(ntm,:); % save the values at those indexes in the cell
end
