Fmuestreo = 200
fout = open("marcasTiempo.txt","w")
with open("marcas.txt") as fin:
	lines = list(map(float, fin.readlines()))
	for muestra in lines:
		fout.write("{}\n".format(muestra/Fmuestreo))

fout.close()

