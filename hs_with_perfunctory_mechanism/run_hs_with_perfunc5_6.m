cp.NP = 100;
cp.D = 1000;
cp.xmin = -100;
cp.xmax = 100;
cp.gen_max = 1000;
cp.max_fes = 3000000;
cp.fl = 0.5;
cp.fs = 0.15;
cp.cr = 0.9;
cp.B = 1.5;
cp.fbias = -450;
cp.HMCR = 0.9;
cp.PAR = 0.3;
cp.bw = 0.01 * 200;

restoredefaultpath
path(path,'../cec2010');
path(path,'../cec2010/datafiles');
fh = fopen('results2010.txt','w');
for func_num = 5:6

   % set the lower and upper bound for each function
   if (func_num == 1 | func_num == 4 | func_num == 7 | func_num == 8 | func_num == 9 | func_num == 12 | func_num == 13 | func_num == 14 | func_num == 17 | func_num == 18 | func_num == 19 | func_num == 20)
      cp.xmin = -100;
      cp.xmax = 100;
   end
   if (func_num == 2 | func_num == 5 | func_num == 10 | func_num == 15)
      cp.xmin = -5;
      cp.xmax = 5;
   end
   if (func_num == 3 | func_num == 6 | func_num == 11 | func_num == 16)
      cp.xmin = -32;
      cp.xmax = 32;
   end

   for run = 1:5
	fitness_val = hs_with_perfunc('benchmark_func',cp,func_num,run);
        fitness_arr(run) = fitness_val;
   end
   fprintf(fh, 'function = %d,min(val) = %f,mean = %d,standard_deviation = %f,\n\n', func_num,min(fitness_arr),mean(fitness_arr),std(fitness_arr));
end

fclose('all');
restoredefaultpath
