function demo()
%Demo Download a json file from Allen Brain Observatory S3 bucket
    strLocalFilename = [tempname, '.json'];
    strURLFilename = 'https://allen-brain-observatory.s3.us-west-2.amazonaws.com/visual-coding-2p/cell_specimens.json';
    C = onCleanup(@(filename) delete(strLocalFilename));
    downloadFile(strLocalFilename, strURLFilename, ShowFilename=true)
end