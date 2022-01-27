

currdir = pwd;

tmp = strsplit(currdir, 'ptutils');

fprintf('\ncalibration functions:\n')
dir([tmp{1} 'ptutils/utils/calibfuncs/*.m']);

