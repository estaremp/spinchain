#Makefile

#Compile

#Define directories
src = src
plot = src/Plotting
bin = bin

#Define variables
main = $(bin)/MAIN.o  
modules = $(src)/modules.f90
plots = $(bin)/*.pyc

#Compilers 
f90comp = gfortran
Pycomp = python2 

#Libraries
libs = -llapack -lblas

#Compile fortran first create bin, then execute
running: $(bin) | exec 

#Execute procedure, first compile programs then plots
exec: $(main) | $(plots)
	$(f90comp) $(libs) $(src)/modules.o $(main) -o run
	@echo "||||Build succesful!||||"

#Create bin
$(bin):
	@echo "|||Creating bin||||"
	mkdir -p $(bin)

#Compile programs, first modules
$(bin)/MAIN.o: $(src)/MAIN.f90 
	@echo "|||Compiling Fortran||||"
	$(f90comp) -c $(modules) -o $(src)/modules.o
	$(f90comp) -c $^ -o $@
	
#Compile python plots
$(bin)/%.pyc: $(plot)/%.py
	@echo "|||Compiling Python||||"
	$(Pycomp) -m py_compile $^ 

#Clean
clean:
	@echo "|||Cleaning|||"
	rm run
	rm -f *.mod
	rm -f $(src)/*.o
#Endof Makefile
