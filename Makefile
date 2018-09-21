#sps
#FC = f77 -g -J
#sun
#FC = f77
#FC est generalement defini par defaut
#FC=ifort
FC = gfortran

FFLAGS = -O

.c.o:
	$(CC) $(CFLAGS) -c $*.c

.f.o:	$(INC)
	$(FC) $(FFLAGS) -c $*.f
	
clean:
	/bin/rm *.o axitra source

axitra :	axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4v3.o reflect5v3.o ff0ad.o source.o fft2cd.o stress_force.o
	$(FC) $(FFLAGS) -o axitra axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4v3.o reflect5v3.o ff0ad.o source.o fft2cd.o stress_force.o
