
include ../make.inc

LINKS = \
-I$(LIB_DIRECTORY)

EXECUTABLE = ./build_lib.exe

type_LagrangeInterpolator.o: \
type_LagrangeInterpolator.f90
	$(FC) $(FFLAGS) -c type_LagrangeInterpolator.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

lagrange_interpolator_library.o: \
type_LagrangeInterpolator.o \
lagrange_interpolator_library.f90
	$(FC) $(FFLAGS) -c lagrange_interpolator_library.f90 $(LINKS) \
	 -o $(OBJ_DIRECTORY)/$@ -J$(LIB_DIRECTORY)

main.o: \
lagrange_interpolator_library.o \
main.f90
	$(FC) $(FFLAGS) -c main.f90 $(LINKS) \
	-o $(OBJ_DIRECTORY)/$@

exec: $(OBJ)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(OBJS) $(LINKS)

lib: 
	$(AR) rc lib$(LIB_NAME).a $(OBJS)
	mv lib$(LIB_NAME).a $(LIB_DIRECTORY)

all:
	$(MAKE) clean
	$(MAKE) exec; $(EXECUTABLE)
	$(MAKE) lib

clean:
	rm -f $(OBJS) $(LIB) $(LIB_DIRECTORY)/*.mod
	rm -f $(EXECUTABLE)