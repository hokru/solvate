
  PROG = Solvate

  OBJS=modules.o solvate.o random.o wpack.o string.o read_options.o molecule.o

  FC = gfortran 
  LINK =  -static -fopenmp
  FLAGS= -O3 -ffree-line-length-none -m64  


#  any blas/lapack will do
#  OPENBLAS=/usr/qc/openblas
#  LIBS= -L$(OPENBLAS) -lopenblas 

  LIBS=-L/usr/lib64 -llapack -lblas


# targets:
.PHONY: all
.PHONY: clean


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FLAGS) -c $< -o $@

$(PROG):$(OBJS) 
	$(FC) $(LINK) $(OBJS) $(LIBS) -o $(PROG)


clean:
	rm -f *.o *.mod $(PROG) 

tar:
	rm -f solv.tgz
	tar -cvzf solv.tgz Makefile *.f90 README.md LICENSE.md database
