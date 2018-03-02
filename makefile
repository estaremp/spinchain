#Makefile

#Compile

#Define directories
src = src
plot = src/Plotting
bin = bin

#Define variables
main = $(bin)/main.o  
const = $(src)/CONSTANTS.f90
param = $(src)/PARAMETERS.f90
depen = $(src)/DEPENDENCIES.f90
plots = $(bin)/*.pyc

#Compilers 
f90comp = gfortran
Pycomp = python2 

#Libraries
libs = -llapack -lblas

#flags
flags = -g -Wall -pedantic -fbounds-check  -Wtabs -fbacktrace -ffpe-trap=invalid,zero,overflow

#Compile fortran first create bin, then execute
running: $(bin) | exec 

#Execute procedure, first compile programs then plots
exec: $(main) | $(plots)
	$(f90comp) $(bin)/c.o $(bin)/p.o $(bin)/d.o $(main) -o run $(libs)
	@echo "***********************"
	@echo "||||Build succesful||||"
	@echo "***********************"


#Create bin
$(bin):
	@echo "***********************"
	@echo "|||Creating bin||||"
	@echo "***********************"
	mkdir -p $(bin)

#Compile programs, first modules
$(bin)/main.o: $(src)/main.f90 
	@echo "***********************"
	@echo "|||Compiling Modules||||"
	@echo "***********************"
	$(f90comp) -c $(const) -o $(bin)/c.o
	$(f90comp) -c $(param) -o $(bin)/p.o
	$(f90comp) -c $(depen) -o $(bin)/d.o
	@echo "***********************"
	@echo "|||Compiling Main||||"
	@echo "***********************"
	$(f90comp) -c $^ -o $@
	
#Compile python plots
$(bin)/%.pyc: $(plot)/%.py
	@echo "***********************"
	@echo "|||Compiling Python||||"
	@echo "***********************"
	$(Pycomp) -m py_compile $^ 

#Clean
clean:
	@echo "|||Cleaning|||"
	rm run
	rm -f *.mod
	rm -f $(src)/*.o
#Endof Makefile
