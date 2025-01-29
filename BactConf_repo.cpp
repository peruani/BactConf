//                                                                                                                                                       
//  To compile: c++ BactConf_repo.cpp -o BactConf_repo -lgsl -lgslcblas -lm -O3	//-static                                                              
//                               
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>       /* srand, rand */
#include <ctime>        /* time */

#define PI 3.14159265358979323846

using namespace std;
using namespace std::chrono;


int main(int argc, const char * argv[]) { //open the main function

//random number generator
    mt19937_64 rng;
    uint64_t timeSeed = chrono::high_resolution_clock::now().time_since_epoch().count(); //it uses time for the seed.
    uniform_real_distribution<double> unif(0, 1); //it generates a random uniform number between zero and one.
    normal_distribution<double> norm(0, 1); //it generates a random number following a normal distribution with average zero and stdv one.
    //Very important that these two lines of the seed are inside the main, otherwise, an error will pop up.
    seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss); //if I want to debug the code, I need to set this ss to a constant value, e.g. 10, so that I get the same set of random numbers and then when I print things it's easier for me to debug.


    //****************************************************************

    int lFlagTraj = 1;

    double v0 = 27.0000000000000; 
    int MaxRW = 10000; //number of walkers
    int MaxEventos = 100000000;
    int MaxSaltos = 30000000.0000000000000; 
    double Dt = 0.001; 
    double D = 9.65; 
    
    double tau0 = exp(1.53);
    double tau_M = 19; //correlation time
    double Delta_n = 1.62;
    double t_CW = 0.4;
    // int max_ContadorTau = 100000; //size array times;
    int save_cada = 10;// Every save_cada integration steps MSD and XYZ (if lFalgTraj=1) is saved  
    int FlagModel = 1;

    double tildeH;
    double H;
    double tildeR;
    double R;
    double omega = 0.0000000000000; 


    if (argc>1) {
	  char * pEnd;   
	  FlagModel 	= atoi(argv[1]); 	// FlagModel=0 ** Poisson **; FlagModel=1 ** BV ** 
	  MaxRW 	= atoi(argv[2]); 	// Number of walkers
	  MaxSaltos 	= atoi(argv[3]);
	  Dt 		= strtod(argv[4], &pEnd);
	  tildeH	= strtod(argv[5], &pEnd); 
	  H = tildeH*v0*2.23;
	  tildeR	= strtod(argv[6], &pEnd);
	  R = tildeR*v0*2.23;
	  if (R<0) {
	  	omega=0;
	  } else {
	  	if (R>0) {
	  		omega = v0/R;
	  	} else {
	  		omega = -1/Dt; 
	  	}
	  }
    } 
    
    int sizePosition = MaxSaltos/save_cada;
    int PosByWalkers = sizePosition*MaxRW;


    //Initialization of everything:
    int i,j; 
    int max_ContadorTau = 100000; //hago 20mil ya que es lo que hace aprox. Eric. //100000;

    double *arr_tau_run = (double*) malloc(max_ContadorTau * sizeof(double));
    double *arr_tau_tum = (double*) malloc(max_ContadorTau * sizeof(double));
    double *arrPositionX = (double*) malloc(sizePosition * sizeof(double));
    double *arrPositionY = (double*) malloc(sizePosition * sizeof(double));
    double *arrPositionZ = (double*) malloc(sizePosition * sizeof(double));
    double *arrPositionSqXPlusSqY = (double*) malloc(sizePosition * sizeof(double));
   double *arrTiempo  = (double*) malloc(sizePosition * sizeof(double));
    

    double *arrTbulk = (double*) malloc(MaxEventos * sizeof(double));
    double *arrTsurface = (double*) malloc(MaxEventos * sizeof(double));


    for (int i = 0; i < max_ContadorTau; i++)
    {
    	arr_tau_run[i]=0;
    	arr_tau_tum[i]=0;
    }

    for (j = 0; j < sizePosition; j++) {
	arrPositionX[j] = 0.0;
	arrPositionY[j] = 0.0;
	arrPositionZ[j] = 0.0;
	arrPositionSqXPlusSqY[j] = 0.0;
	arrTiempo[j] = 0.0;
    }
    

//*****************************Output files**************************************    
    
    // ************ File names ***********************
    char filename[100];
    strcpy(filename,"_Model_"); strcat(filename,argv[1]);
    strcat(filename,"_RW_"); strcat(filename,argv[2]);
    strcat(filename,"_T_"); strcat(filename,argv[3]);
    strcat(filename,"_Dt_"); strcat(filename,argv[4]);
    strcat(filename,"_tH_"); strcat(filename,argv[5]);
    strcat(filename,"_tR_"); strcat(filename,argv[6]);
    strcat(filename,".dat");
  
 
    char filename_XYZ[100];
    strcpy(filename_XYZ,"XYZ"); strcat(filename_XYZ,filename);
    ofstream toXYZ(filename_XYZ);          			 // open file 
    toXYZ.precision(8);
    toXYZ.setf(ios::scientific,ios::floatfield);

    char filename_MSD[100];
    strcpy(filename_MSD,"MSD"); strcat(filename_MSD,filename);
    ofstream toMSD(filename_MSD);          			  
    toMSD.precision(8);
    toMSD.setf(ios::scientific,ios::floatfield);

    double deltaXp = 0;
    double p_tumb;
    double tt;
    int n_CW;
    
    double a = 1 - Dt/tau_M; 
    double b = sqrt(Dt)*sqrt(2/tau_M); 
    double inv_Dt_tau0 = Dt/tau0;

    int save_tbulk = 0;
    int save_tsurf = 0;
   
    int RW;
    double e1; 
    double e2; 
    double e3; 
    double e1_new;
    double e2_new;
    double e3_new;
    double norm_e;
    int coeffOmega=1;

    double eta1, eta2, eta3;
    double D_theta = 0.0;


    int contadorTau = 0;
    double accumulate = 0;
    double orientCorr;
    double dir1_before, dir2_before, dir3_before;
    double dir1_after, dir2_after, dir3_after;
    double inner_product;
    int cont_tumbling = 0;
    int contadorRun;
    int contadorRunBulk;
    int contadorRunSurf;
    int contadorTum;
    int maxCont_tumbling = 10000;
    int mover_tiempo;
    double theta,alpha;
    double X,Y,Z,X_ini,Y_ini,Z_ini;
    double tbulk, tsurface;
    int lbulk_event, lsurface_event;
    int time_cont, almaceno;
    int lFlagSup = 1;
    int lFlagSigoEnSup = 1;

    int toco_abajo;
    int toco_arriba;
    double sum_touchUpDown;
    double av_SumtouchUpDown=0.0;
    double ratio;
    double av_ratio=0.0;
    double var_ratio=0.0;
    double std_ratio;
    double difference;
    double av_difference=0.0;
    double var_difference=0.0;
    double std_difference;
    

    const int num_bins = 10000; 
    double bin_size = H/num_bins;
    int my_int;
    double eje_z[num_bins+1];
    for (i=0; i<num_bins+1; i++) eje_z[i] = i*bin_size;
    double PdeZ[num_bins+1];
    for (i=0; i<num_bins+1; i++) PdeZ[i] = 0.0;

    double arrInnerProduct[maxCont_tumbling];
    for (i=0; i<maxCont_tumbling; i++) arrInnerProduct[i] = 0.0;
    
    contadorRun=0;
    contadorRunBulk=0;
    contadorRunSurf=0;

    contadorTum=0;
    int clock=0;
    int clock2=0;

    //===================== Loop over walkers 
    for (RW = 0; RW < MaxRW; RW++) { 
           
	    theta = unif(rng)*2*PI;
	    e1 = cos(theta);
	    e2 = sin(theta);
	    e3 = 0;
	    X = 0.0;
	    X_ini = 0.0;
	    Y = 0.0;
	    Y_ini = 0.0;
	    Z = 0.0; 
	    Z_ini = 0.0; 
	   
	    if ((Z<=0) || (Z>=H)) {    	
	  	lFlagSup = 1; 
	    } else {
	    	lFlagSup = 0;
	    }
	    
	    deltaXp = 0;
	
	    clock=0;
	    clock2=0;

	    time_cont = 0;
	    almaceno = 0;
	    
	
    	    //===================== Loop over time 
	    while (time_cont<MaxSaltos) { 

		if (lFlagSup == 0) { 
			tbulk = tbulk + Dt;
			X = X + v0*e1*Dt;
			Y = Y + v0*e2*Dt;
			Z = Z + v0*e3*Dt;

	        	if (Z <= 0) {
	            		lFlagSup = 1;    
	            		Z = 0;
			    	e3 = 0;
			    	norm_e = sqrt(e1*e1 + e2*e2 + e3*e3);

				if (norm_e <= 0.00001) {
					theta = unif(rng)*2*PI;
		    			e1 = cos(theta);
		    			e2 = sin(theta);
					norm_e = 1;
				} 
				else {
					e1 = e1/norm_e;
					e2 = e2/norm_e;
					theta=atan2(e2,e1);
				}    
		        } 
		        if (Z >= H) { 
		            lFlagSup = 1;
			    Z = H;
			    e3 = 0;
			    norm_e = sqrt(e1*e1 + e2*e2 + e3*e3);
                            if (norm_e <= 0.00001) {
				    theta = unif(rng)*2*PI;
            			    e1 = cos(theta);
            			    e2 = sin(theta);
            		     }
            		     else {
				    e1 = e1/norm_e;
			            e2 = e2/norm_e;
				    theta=atan2(e2,e1);
                            }

		        }
		        		     
		}
		else {
                        tsurface = tsurface + Dt;
                        if (omega>=0) {
                        	X = X + v0*cos(theta)*Dt;
                        	Y = Y + v0*sin(theta)*Dt;
                        }

                        if (Z>=H) {
                              coeffOmega=-1; 
                        }
                        if (Z<=0) {
                                coeffOmega=1;
                        }
                        theta = theta + Dt*coeffOmega*omega;
               } 
		mover_tiempo=1;

		clock++;
		clock2++;
		
		if (FlagModel>0) {
			deltaXp = a*deltaXp + b*norm(rng); 
			p_tumb = exp(Delta_n*deltaXp)*inv_Dt_tau0;
		} else {		
			p_tumb = Dt/6.74;
		}

		if (unif(rng) < p_tumb) { //open if I tumble
			if (contadorRun < max_ContadorTau) { 				
				contadorRun++;
				arr_tau_run[contadorRun]=clock*Dt;

			}
			clock=0;

		        tt = (-t_CW)*log(1-unif(rng)); 
		        n_CW = max(1.0, floor(tt/Dt)); //Time it stays in CW mode (CW = tumble mode)
		     

		        norm_e = sqrt(e1*e1 + e2*e2 + e3*e3);
			e1 = e1/norm_e;
			e2 = e2/norm_e;
			e3 = e3/norm_e; 
			e1_new = e1;
			e2_new = e2;
			e3_new = e3;	
		        

			dir1_before = e1;
			dir2_before = e2;
			dir3_before = e3;
			
			//=== reorientation of e versor
			for (j = 0; j < n_CW; j++) { 
				deltaXp = a*deltaXp + b*norm(rng);
				eta1 = norm(rng); eta2 = norm(rng); eta3 = norm(rng);
				e1_new = e1 + Dt*(-2)*D*e1 + sqrt(2*D*Dt)*((1 - e1*e1)*eta1 - e1*e2*eta2 - e1*e3*eta3);
				e2_new = e2 + Dt*(-2)*D*e2 + sqrt(2*D*Dt)*((1 - e2*e2)*eta2 - e2*e3*eta3 - e2*e1*eta1);
				e3_new = e3 + Dt*(-2)*D*e3 + sqrt(2*D*Dt)*((1 - e3*e3)*eta3 - e3*e1*eta1 - e3*e2*eta2);
				norm_e = sqrt(e1_new*e1_new + e2_new*e2_new + e3_new*e3_new);
				e1 = e1_new/norm_e;
				e2 = e2_new/norm_e;
				e3 = e3_new/norm_e;
		        } 
			

			if (lFlagSup == 1) { 
				if (Z <= 0) { 
					// Bottom surface
					if (e3>0) {						
						save_tbulk++;
						lFlagSup = 0;
					}
				}
				else {
					// Upper surface
					if (e3<0) {						
						save_tbulk++;
						lFlagSup = 0;
					}
				}
				if (lFlagSup == 1) {
					// if lFlagSup has NOTchanged, then RW remains on the surface
					e3 = 0;
			    		norm_e = sqrt(e1*e1 + e2*e2 + e3*e3);
					if (norm_e <= 0.00001) {
						theta = unif(rng)*2*PI;
	    					e1 = cos(theta);
	    					e2 = sin(theta);
						norm_e = 1;
					} 
					else {
						e1 = e1/norm_e;
						e2 = e2/norm_e;
						theta=atan2(e2,e1);
					}
				}    		
				
			} //close if surface
			

			dir1_after = e1;
			dir2_after = e2;
			dir3_after = e3;
			inner_product = dir1_before*dir1_after + dir2_before*dir2_after + dir3_before*dir3_after;
			accumulate = accumulate + inner_product;

			if (cont_tumbling < maxCont_tumbling) {
				arrInnerProduct[cont_tumbling] = inner_product;
			}
			cont_tumbling++;
			

			
			//=== move tiempo + n_CW !!!!
			mover_tiempo=mover_tiempo+n_CW;

		} // close tumbling event 

		for (int kk=0;kk<mover_tiempo;kk++) {
			time_cont++;
			if ((time_cont % save_cada == 0) && (almaceno < sizePosition)) { 

		        	arrPositionSqXPlusSqY[almaceno] = arrPositionSqXPlusSqY[almaceno] + (X - X_ini) * (X - X_ini) + (Y - Y_ini) * (Y - Y_ini);
			
				if (lFlagTraj==1) {
					arrPositionX[almaceno] = arrPositionX[almaceno] + X;
					arrPositionY[almaceno] = arrPositionY[almaceno] + Y;
					arrPositionZ[almaceno] = arrPositionZ[almaceno] + Z;
				}
			
		        	arrTiempo[almaceno] = time_cont * Dt; //tbulk + tsurface
		        	almaceno++;
			}
		}

	  } //close while run times OR close while eventos, según lo que estemos haciendo.
	  	  
	} //close for walkers loop
	

	
	
  
    for (j = 0; j < sizePosition; j++) { //open average over positions' loop

	arrPositionSqXPlusSqY[j] = arrPositionSqXPlusSqY[j] / MaxRW;
	toMSD << arrTiempo[j] << "\t"<< arrPositionSqXPlusSqY[j] << endl;

	if (lFlagTraj==1) {
		arrPositionX[j] = arrPositionX[j] / MaxRW;
		arrPositionY[j] = arrPositionY[j] / MaxRW;
		arrPositionZ[j] = arrPositionZ[j] / MaxRW;
		toXYZ << arrTiempo[j] << "\t"<< arrPositionX[j] << "\t" << arrPositionY[j] << "\t" << arrPositionZ[j] << endl;
	}

    } //close average over positions' loop



    

return 0;
} //close the main function



































