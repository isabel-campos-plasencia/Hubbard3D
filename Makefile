CFLAGS= -g -Wall                             #para debug  
CFLAGS=  -O3 -funroll-loops -finline-functions -pg #optimizacion maxima

OBJ   = hm3d.o hm3d_matrix.o hm3d_ini.o hm3d_stable.o hm3d_green.o hm3d_update.o hm3d_measure.o nrutil.o hm3d_nr.o 

iso: $(OBJ)
	gcc $(CFLAGS) $(OBJ) -lm -o Hubbard_3d
.c.o:
	gcc $(DEFINES) -c $(CFLAGS) $<

clean: 
	/bin/rm -f $(OBJ) 

#		*Individual File Dependencies*

hm3d.o: hm3d.c hm3d.h 

hm3d_matrix.o: hm3d_matrix.c hm3d.h

hm3d_ini.o: hm3d_ini.c  hm3d.h 

hm3d_stable.o: hm3d_stable.c hm3d.h

hm3d_green.o: hm3d_green.c hm3d.h

hm3d_update.o: hm3d_update.c hm3d.h

hm3d_measure.o: hm3d_measure.c hm3d.h

nrutil.o: nrutil.c hm3d.h

hm3d_nr.o: hm3d_nr.c hm3d.h



