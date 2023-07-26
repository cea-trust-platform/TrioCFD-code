# Sensitivity
########## Sensitivity section ##########
def figures(variable, cross_sec, direction):
	test_case = "Sensitivity_analysis_for_Steady_Obstacle_"
	confidence = 0.95 #level of confidence for the interval
	std_dev = 0.0075 #std deviation of the input

	# MC
	N = 1204 # number of MC simulations
	test_case_MC = "Obstacle_"

	import numpy as np
	import matplotlib.pyplot as plt
	x = []
	y = []
	with open(test_case+variable+'_'+cross_sec+'.son', 'r') as f:
		for line in f:
			state=[]
			if line[0] != '#':
				for num in line.split(' '):
					state.append(float(num))
			else:
				text_old = " "
				for text in line.split(' '):
					if text_old == "x=":
						x.append(float(text))
					if text_old == "y=":
						y.append(float(text))
					text_old = text

	state = np.delete(state, 0, 0)
	with open(test_case+variable+'_'+cross_sec+'_SENS'+'.son', 'r') as f:
		for line in f:
			sens = []
			if line[0] != '#':
				for num in line.split(' '):
					sens.append(float(num))
	sens = np.delete(sens, 0)

	if np.var(y) < np.var(x):
		abscissa = x
		xlabel = "x"
		section = "y"
		section_lev = str(round(np.mean(y), 2))
	else:
		abscissa = y
		xlabel = "y"
		section = "x"
		section_lev = str(round(np.mean(x), 2))

	variance = np.abs(sens)*std_dev
	upper_bound = state+1/np.sqrt(1-confidence)*variance
	lower_bound = state-1/np.sqrt(1-confidence)*variance
	if variable == "P":
		average = state
		var = "the pressure"

	if variable == "U":
		if direction == "x":
			average = state[0::2]
			lower_bound = lower_bound[0::2]
			upper_bound = upper_bound[0::2]
			variance = variance[:][0::2]
			var = "$u_x$"
		if direction == "y":
			average = state[1::2]
			lower_bound = lower_bound[1::2]
			upper_bound = upper_bound[1::2]
			variance = variance[:][1::2]
			var = "$u_y$"


	########## MC section ##########

	average_MC = np.loadtxt('../Resultats_MC/average_MC_'+variable+'_'+cross_sec+'.txt', unpack=True)
	variance_MC = np.loadtxt('../Resultats_MC/variance_MC_'+variable+'_'+cross_sec+'.txt', unpack=True)  

	if variable == "U":
		if direction == "x":
			average_MC = average_MC[:][0::2]
			variance_MC = variance_MC[:][0::2]
			var = "$u_x$"
			var_for_fig = 'Ux'
		if direction == "y":
			average_MC = average_MC[:][1::2]
			variance_MC = variance_MC[:][1::2]
			var = "$u_y$"
			var_for_fig = 'Uy'
	if variable == "P":
		var_for_fig='P'
	########## Plots ##########



	fig, ax = plt.subplots()
	ax.plot(abscissa, average)
	ax.plot(abscissa, average_MC,   c = 'red')
	ax.legend(['Average SA','Average MC'])
	plt.savefig('../Resultats_MC/Average_'+var_for_fig+'_'+cross_sec+'.png')

	fig, ax = plt.subplots()
	ax.plot(abscissa, variance)
	ax.plot(abscissa, np.sqrt(variance_MC),  c = 'red')
	ax.legend([' Standard deviation SA',' Standard deviation MC'])
	plt.savefig('../Resultats_MC/Variance_'+var_for_fig+'_'+cross_sec+'.png')
	plt.close ()

	fig, ax = plt.subplots()
	ax.plot(abscissa, upper_bound,  c = 'red')
	ax.plot(abscissa, average)
	ax.plot(abscissa, lower_bound,  c = 'red')
	ax.legend(['Upper bound', 'Average SA','Lower bound'])
	plt.savefig('../Resultats_MC/CI_'+var_for_fig+'_'+cross_sec+'.png')
	plt.close ()
figures("U", "H2", "x")
figures("U", "H2", "y")

figures("U", "H3", "x")
figures("U", "H3", "y")

figures("U", "V6", "x")
figures("U", "V6", "y")

figures("U", "V1", "x")
figures("U", "V1", "y")

figures("P", "H2", "x")
figures("P", "H3", "x")
figures("P", "V6", "x")
figures("P", "V1", "x")

