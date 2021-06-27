all: main.o hcubature.o  
	cc  main.o hcubature.o
main.o: main.c
	cc main.c
hcubature.o: hcubature.c
	cc hcubature.c 
