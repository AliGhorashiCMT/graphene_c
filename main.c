#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "cubature.h"

double pi = 3.14159265;
double a = 2.46;
//double a1[2] = {a, 0};
//double a2[2] = {-a/2, a*sqrt(3)/2};
//double b1[2] = {2*pi/a, 2*pi/(sqrt(3)*a)};
//double b2[2] = {0, 4*pi/(a*sqrt(3))};

int plasmon(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    double q = *((double *) fdata); // we can pass σ via fdata argument
    fdata = fdata + sizeof(double);
    double omega =  *((double *) fdata);
    // printf("q = %f\n", q );
    // printf("omega = %f\n", omega);
    double k = x[0];
    double theta = x[1];
    if(6*k > 1){
        fval[0] = 0;
        return 0;
    }
    else{
        //fval[0] = k;
         double ek = 6*k;
         double ekq = 6*sqrt(k*k+q*q+2*k*q*cos(theta));
         fval[0] = -(k/(pi*pi))*creal(2*(ekq-ek)/((ekq-ek)*(ekq-ek)-(omega+.001*I)*(omega+.001*I)));
    }
    // compute the output value: note that fdim should == 1 from below
    return 0; // success*
}

double plasmon_calculator(double q, double omega){
    double val;
    double err;
    double qomega[] ={q, omega};
    void * pqomega = qomega;
    double xmin[2] = {0, 0}; double xmax[2] = {1, 2*pi};
    hcubature(1, plasmon, pqomega, 2, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
    double epsilon = 1-(90.5/qomega[0])*val;
    //printf("Epsilon = %f\n", epsilon);
    //printf("Computed integral = %0.10g +/- %g\n", val, err);
    return epsilon;
}

double in_unit_cell(double x, double y){
	double d0 = pow(x, 2)+pow(y,2);
	double d1 = pow(x-a, 2)+pow(y, 2);
	double d2 = pow(x+a, 2)+pow(y, 2);
	double d3 = pow(x-a/2.00, 2)+pow(y-a*sqrt(3)/2.00, 2);
        double d4 = pow(x+a/2.00, 2)+pow(y+a*sqrt(3)/2.00, 2);
        double d5 = pow(x-a/2.00, 2)+pow(y+a*sqrt(3)/2.00, 2);
        double d6 = pow(x+a/2.00, 2)+pow(y-a*sqrt(3)/2.00, 2);
	if ( (d0 < d1 ) && (d0 < d2) && (d0 < d3) && (d0 < d4) && (d0 < d5) && (d0 < d6) ){
		return 1;
	}
	else{
		return 0;
	}
}

double in_brillouin_zone(double kx, double ky){
	double b = 4*pi/(a*sqrt(3));
	double d0 = pow(kx, 2) + pow(ky, 2);
	double d1 = pow(kx, 2) + pow(ky-b, 2);
	double d2 = pow(kx, 2) + pow(ky+b, 2);
	double d3 = pow(kx-b*sqrt(3)/2.00, 2) + pow(ky-b/2.00, 2);
	double d4 = pow(kx+b*sqrt(3)/2.00, 2) + pow(ky+b/2.00, 2);
	double d5 = pow(kx-b*sqrt(3)/2.00, 2) + pow(ky+b/2.00, 2);
	double d6 = pow(kx+b*sqrt(3)/2.00, 2) + pow(ky-b/2.00, 2);

        if ( (d0 < d1 ) && (d0 < d2) && (d0 < d3) && (d0 < d4) && (d0 < d5) && (d0 < d6) ){
                return 1;
        }
        else{
                return 0;
        }
}

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
	double meshu = atof(argv[3]);
	double meshb = atof(argv[4]);
	double meshplas = atof(argv[5]);

	dos2(mesh, width);
	
	FILE *fpcell;
        fpcell = fopen("./unitcell.txt", "w"); 
	FILE *fpbcell;
        fpbcell = fopen("./brillouinzone.txt", "w");

	for(int i=1; i<meshu; i++){
		double x = 5.00*(double)(i-meshu/2.000)/meshu;
		for(int j=1; j< meshu; j++){
                	double y = 5.00*(double)(j-meshu/2.000)/meshu;
			fprintf(fpcell, "%f\t", in_unit_cell(x, y));
		}
		fprintf(fpcell, "\n");
	}
	fclose(fpcell);

        for(int i=1; i<meshb; i++){
                double kx = 2*K[0]*(double)(i-meshb/2.00)/meshb;
                for(int j=1; j<meshb; j++){
                        double ky = 2*K[0]*(double)(j-meshb/2.00)/meshb;
                        fprintf(fpbcell, "%f\t", in_brillouin_zone(kx, ky));
                }
                fprintf(fpbcell, "\n");
        }
        fclose(fpbcell);
        FILE *epsilonfile;
        epsilonfile = fopen("./epsilon.txt", "w");


        for(int i=1; i<meshplas; i++){
            printf("Iteration number: %d\n of %f ", i, meshplas );
            double q = ((double) i)/(meshplas)*(2.00/6.000);
	    for(int j=1; j<meshplas; j++){
                 double omega = ((double) j)/(meshplas)*(3.00);
		  double epsilon = plasmon_calculator(q, omega);
	         fprintf(epsilonfile, "%f\t", epsilon);
            }
            fprintf(epsilonfile, "\n");
        }
        fclose(epsilonfile);
	return 0;
}
