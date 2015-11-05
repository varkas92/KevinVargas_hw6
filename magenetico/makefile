kinetic_energy = 0.00001
alpha = 30

particle_xy.pdf: trayectoria_$(kinetic_energy)_$(alpha).dat particle_in_field.py
	python particle_in_field.py $(kinetic_energy) $(alpha)

trayectoria_$(kinetic_energy)_$(alpha).dat: particle_in_field.x
	./particle_in_field.x $(kinetic_energy) $(alpha) > trayectoria_$(kinetic_energy)_$(alpha).dat

particle_in_field.x: particle_in_field.c
	cc particle_in_field.c -lm -o particle_in_field.x
