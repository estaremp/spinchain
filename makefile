#Makefile

#Compile

#Define directories
src = src
plot = src/Plotting
bin = bin

#Define variables
program = $(bin)/main.o
dependencies = $(src)/modules.f90
plots = $(bin)/*.pyc

#Compilers 
f90comp = gfortran
Pycomp = python2 

#Libraries
libs = -llapack -lblas

#Compile fortran first create bin, then execute
running: $(bin) | exec 

#Execute procedure, first compile programs then plots
exec: $(program) | $(plots)
	$(f90comp) $(libs) $(src)/modules.o $(program) -o run
	@echo "||||Build succesful!||||"

#Create bin
$(bin):
	@echo "|||Creating bin||||"
	mkdir -p $(bin)

#Compile programs, first modules
$(bin)/%.o: $(src)/%.f90 
	@echo "|||Compiling Fortran||||"
	$(f90comp) -c $(dependencies) -o $(src)/modules.o
	$(f90comp) -c $^ -o $@
	
#Compile python plots
$(bin)/%.pyc: $(plot)/%.py
	@echo "|||Compiling Python||||"
	$(Pycomp) -m py_compile $^ 

#Clean
clean:
	@echo "|||Cleaning|||"
	rm run
	rm -f $(program)
	rm -f $(src)/modules.o
	rm -f *.mod
#Endof Makefile
