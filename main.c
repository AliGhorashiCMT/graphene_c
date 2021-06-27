#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>


double pi = 3.14159265;
double a = 2.46;
//double a1[2] = {a, 0};
//double a2[2] = {-a/2, a*sqrt(3)/2};
//double b1[2] = {2*pi/a, 2*pi/(sqrt(3)*a)};
//double b2[2] = {0, 4*pi/(a*sqrt(3))};


double full_dispersion(double kx, double ky){
	double a_nn = 1.42; //Angstroms
	double hopping = 2.8; //eV
	double complex add_hoppings = cexp(I*(a_nn)*(ky))+cexp(I*(a_nn)*(-ky/2.000+kx*sqrt(3.000)/2.000))+cexp(a_nn*I*(-ky/2.000-kx*sqrt(3.000)/2.00));
	return hopping*csqrt(creal(add_hoppings)*creal(add_hoppings)+cimag(add_hoppings)*cimag(add_hoppings));
}
int dos2(double mesh, double width){
	printf("calculating dos\n");
        FILE *fpdos;
        fpdos = fopen("./dos.txt", "w");
        int idx;
	printf("Using mesh %f and histogram width %f\n", mesh, width);
	
        double b1[2] = {2*pi/a, 2*pi/(sqrt(3)*a)};
        double b2[2] = {0, 4*pi/(a*sqrt(3))};
	int larray;
	larray = round(10*width);
	double dosarray[2000] = {0.00};        
	//float dmesh = (float) mesh;
	for (int i=1; i<mesh; i++){
             //printf("%d", i);   
	     for (int j=1; j<mesh; j++){
                        double kx = b1[0]*((double)i)/((double) mesh)+b2[0]*((double)j)/((double) mesh);
                        double ky = b1[1]*((double)i)/((double) mesh)+b2[1]*((double)j)/((double) mesh);
                        double energy = full_dispersion(kx, ky);
                        idx = round(energy*width);
                        dosarray[idx+1000] += width/(mesh*mesh); //(width)*1.0/(mesh*mesh);
                        dosarray[1000-idx] += width/(mesh*mesh); //(width)*1.0/(mesh*mesh);
		}
        }	
	double sum=0;
        for(int i=0; i<2000; i++){
                fprintf(fpdos, "%f\n", dosarray[i]);
        	sum += 1/width*dosarray[i];
	}
	printf("Normalization of density of states: %f", sum);
	fclose(fpdos);
return 0;
}

int main(int argc, char * argv[]){
        double a1[2] = {a, 0};
        double a2[2] = {-a/2, a*sqrt(3)/2};
        double b1[2] = {2*pi/a, 2*pi/(sqrt(3)*a)};
        double b2[2] = {0, 4*pi/(a*sqrt(3))};

	FILE *fp; 
	// K = 2/3*b1-1/3b2
	double K[2]={0, 0};
	double M[2]={0, 0};
	// printf("%p", K);
	/*
	for(int i =0; i< 2; i++){
		printf("%f\t", *(K+i));
	}
	*/
	printf("b1= %f\n\n",2.0/3.0*b1[0]);
	for(int i=0; i<2; i++){
		printf("i=%d", i);
		printf("%f\n", 2.0/3.0*b1[i]-1.0/3.0*b2[i]);
		K[i] += 2.0/3.0*b1[i]-1.0/3.0*b2[i];
		M[i] += 1.0/2.0*b1[i];
	}
	for(int i =0; i< 2; i++){
                printf("K[%d] = %f\t \n", i, *(K+i));
        }
	fp = fopen("./dispersion.txt", "w");
	for(int i=0;i<100; i++){
		double kx = ((double)i)/(100.0)*M[0];
		double ky = ((double)i)/(100.0)*M[1];
		fprintf(fp, "%f\t%f\n", -full_dispersion(kx, ky), full_dispersion(kx, ky));
	}
        for(int i=0;i<100; i++){
                double kx = M[0]+((double)i)/(100.0)*(K[0]-M[0]);
                double ky = M[1]+((double)i)/(100.0)*(K[1]-M[1]);
                fprintf(fp, "%f\t%f\n", -full_dispersion(kx, ky), full_dispersion(kx, ky));
        }
        for(int i=0;i<100; i++){
                double kx = K[0]-((double)i)/(100.0)*K[0];
                double ky = K[1]-((double)i)/(100.0)*K[1];
                fprintf(fp, "%f\t%f\n", -full_dispersion(kx, ky), full_dispersion(kx, ky));
        }
	fclose(fp);
	//printf("%p", K);
	printf("E(Gamma) = %f\n", full_dispersion(0, 0));
	printf("E(K)= %f\n", full_dispersion(K[0], K[1]));
	printf("E(M)= %f\n", full_dispersion(M[0], M[1]));
	
	double mesh = atof(argv[1]);
	double width = atof(argv[2]);
	dos2(mesh, width);

	return 0;
}
