% uncomment a location below -- adapt to local installation

 location = 'home_matteo';
% location = 'home_luca';
% location = 'work_matteo';

restoredefaultpath

if strmatch(location,'home_matteo','exact')    
    dir1='C:\E\dynare\4.3.3\matlab';
    dir2='C:\E\occbin\occbin_20140630\toolkit_files';
    dir3='C:\E\occbin\toolkit_files_private';
    rmpath 'C:\Program Files (x86)\MATLAB\R2009b\toolbox\ident\idobsolete'
elseif strmatch(location,'home_luca','exact')
    dir1='/Applications/dynare/4.3.3/matlab';
    dir2='../occbin_20140630/toolkit_files';
    dir3='../occbin_20140630/toolkit_files_private';
elseif strmatch(location,'work_matteo','exact')  
    dir1='C:\E\dynare\4.3.1\matlab';
    dir2='C:\E\occbin\occbin_20140630\toolkit_files';
    dir3='C:\E\occbin\toolkit_files_private';
    
else 
    error('Specify path to Dynare installation')
end

path(dir1,path);
path(dir2,path);
path(dir3,path);


dynare_config

