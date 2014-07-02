cp.NP = 100;
cp.D = 100;
cp.xmin = -100;
cp.xmax = 100;
cp.gen_max = 1000;
cp.fl = 0.5;
cp.fs = 0.15;
cp.cr = 0.9;
cp.B = 1.5;
cp.fbias = -450;

restoredefaultpath
path(path,'cec2005');
adecm('benchmark_func',cp,1)
restoredefaultpath
