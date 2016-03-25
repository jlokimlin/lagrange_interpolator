
include ../make.inc

LINKS = \
-L$(LIB_DIRECTORY) -llagrange_interpolator \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./test.exe

test.o: \
test.f90
	$(FC) $(FFLAGS) -c test.f90 $(LINKS)\
	 -o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

run:
	make clean
	make exec; $(EXECUTABLE)
	 
clean:
	rm -f $(OBJS) $(EXECUTABLE)