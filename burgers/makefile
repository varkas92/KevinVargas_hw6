burgers.gif: salida.dat burgers.py
	python burgers.py

salida.dat: burgers
	./burgers > salida.dat

burgers: burgers.c
	cc burgers.c -o burgers -lm
