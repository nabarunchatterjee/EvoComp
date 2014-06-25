
function fitness = fitness_evaluation(vect,fnum)

	a1 = vect(1);
	w1 = vect(2);
	a2 = vect(3);
	w2 = vect(4);
	a3 = vect(5);
	w3 = vect(6);
	theta = 2 * pi / 100;
	fitness = 0;
	for t = 0:100
		fitness_vect = a1 * sin(w1 * t * theta + a2 * sin( w2 * t * theta + a3 * sin(w3 * t * theta)));
		fitness_targ = sin(5 * t * theta - 1.5 * sin(4.8 * t * theta + 2 * sin(4.9 * t * theta)));
		fitness = fitness + (fitness_vect - fitness_targ) * (fitness_vect - fitness_targ);
	end
end
