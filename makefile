### Nombre del compilador
ICC=icc 
#ICC=icpc # descomenta esta opcion si tienes Intel C++ compiler 

### Opciones del compildor
OPCS=-g -c -O3
#OPCS=-g -c -O3 -ipo # descomenta esta opcion si tienes Intel C++ compiler 

### Objetos a compilar
#OBJS= main.o mensajes.o distr.o ftemps.o media.o colision.o termo.o libre.o salida.o fteor.o aleat.o
OBJS= main.o SetupJob.o SingleStep.o termo.o salida.o aleat.o

LJMD: $(OBJS)
	$(ICC) $(OBJS) -o LJMD 

main.o: main.c params.h def_functions.h
	$(ICC) $(OPCS) main.c

SetupJob.o: SetupJob.c params.h def_functions.h
	$(ICC) $(OPCS) SetupJob.c

SingleStep.o: SingleStep.c params.h def_functions.h
	$(ICC) $(OPCS) SingleStep.c

termo.o: termo.c params.h def_functions.h
	$(ICC) $(OPCS) termo.c

salida.o: salida.c params.h def_functions.h
	$(ICC) $(OPCS) salida.c

aleat.o: aleat.c params.h def_functions.h
	$(ICC) $(OPCS) aleat.c

all: LJMD $(OBJS)

.PHONY: clean
clean: 
	rm LJMD $(OBJS)

plot: 
	gnuplot -persist fig.plt

compare:
	gnuplot -persist fig2.plt
