#ifndef Models_v6_CPP
#define Models_v6_CPP

#include <vector>
#include <cmath>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <stdlib.h>
#include <algorithm> 
/*
//************************
// A320 Fuel Consumption
//************************
double C_f1 = 0.63330E+00;   //[Kg/(min*kN)]
double C_f2 = 0.85903E+03;   //[knots]
double C_f3 = 0.91340E+01;   //[kg/min]
double C_f4 = 0.79668E+05;   //[feet]
double C_fcr = 0.95423E+00;  //[dimentionless]
//************************
// A320 Engine Thrust
//************************
double C_T_c1 = 0.14104E+06;       //[Newton]
double C_T_c2 = 0.48917E+05;       //[feet]
double C_T_c3 = 0.65004E-10;       //[1/feet^2]
double C_T_c4 = 0.99797E+01;       //[K]
double C_T_c5 = 0.80268E-02;       //[1/K]
double C_T_des_low = 0.27207E-01;  //[dimentionless]
double C_T_des_high = 0.45711E-01; //[dimentionless]
double Hp_des = 0.12398E+05;       //[feet]
double C_T_des_app = 0.13981E+00;  //[dimentionless]
double C_T_des_ld = 0.34750E+00;   //[dimentionless]
double V_des_ref = 0.31000E+03;    //[knots]
double M_des_ref = 0.78000E+00;    //[dimentionless]
*/

//************************************
// A320 Fuel Consumption Coefficients
//************************************
const double C_f1 = (0.63330E+00)*(1.0/60.0);     //[Kg/(min*kN)]-->(1/60) sec/km
const double C_f2 = (0.85903E+03)*(0.514444E-03);   //[knots]   -->   --> Km/s
const double C_f3 = (0.91340E-02)*(1.0/60.0);     //[kg/min]  --> *(1/60)*10^(-3)  --> tons/sec
const double C_f4 = (0.79668E+05)*(0.0003048);   //[feet]    --> *(1/3280.8)     --> Km
const double C_fcr = (0.95423E+00);               //[dimentionless]
//*********************************
// A320 Engine Thrust Coefficients
//*********************************
const double C_T_c1 = 0.14104E+00;                        //[Newton]   -->*10^(-6) -->(ton*Km)/s^2   
const double C_T_c2 = (0.48917E+05)*(0.0003048);         //[feet]     --> *(1/3280.8)  --> Km
const double C_T_c3 = (0.65004E-10)/(0.0003048*0.0003048);//[1/feet^2] --> *(1/3280.8*3280.8) -->1/(Km)^2
const double C_T_c4 = 0.99797E+01;                        //[K]
const double C_T_c5 = 0.80268E-02;                        //[1/K]
const double C_T_des_low = 0.27207E-01;                   //[dimentionless]
const double C_T_des_high = 0.45711E-01;                  //[dimentionless]
const double H_p_des = (0.12398E+05)*(0.0003048);         //[feet]    --> *(1/3280.8)     --> Km
const double C_T_des_app = 0.13981E+00;                   //[dimentionless]
const double C_T_des_ld = 0.34750E+00;                    //[dimentionless]
const double V_des_ref = (0.31000E+03)*(0.514444E-03);      //[knots]   --> Km/s
const double M_des_ref = 0.78000E+00;                     //[dimentionless]
const double C_T_CR = 0.95;
//*********************************
// Atmosphere Model Parameters
//*********************************
//const double To = 288.15;                       //[K] 
//const double Po = 101325.0;                     //[Pa]              
//const double g = 9.80665;                      //[m/s^2]            *(10^-3)  -->[km/s^2] 
//const double beta_T = -0.0065;                  //[k/m]             *(10^3)   -->[k/km]
//const double R = 287.05287;                     //[m^2/(k*s^2)]     *(10^-6)  -->[km^2/(k*s^2)]
//const double H_p_trop = 11.0;                   //[km]
//const double T_ISA_trop = 216.65;               //[K]
//const double P_trop = 22632.0;                   //[Pa] 
//const double k = 1.4;
//*********************************
// Atmosphere Model Parameters
//*********************************
const double To = 288.15;                       //[K] 
const double Po = 101325.0;                     //[Pa]              
const double g = 9.80665E-03;                   //[km/s^2] 
const double beta_T = -0.0065E03;               //[k/km]
const double R = 287.05287E-06;                 //[km^2/(k*s^2)]
const double H_p_trop = 11.0;                   //[km]
const double T_ISA_trop = 216.65;               //[K]
const double P_trop = 22632.0;                   //[Pa] 
const double k = 1.4;
//********************************************
//       Fuel Consumption [Jet Only]
//********************************************

double n(double V_TAS)
{
	return C_f1 * (1.0 + (V_TAS/C_f2));
}

double f_nom(double V_TAS, double Thrust)
{
	return n(V_TAS) * Thrust;
}

double f_min(double H_p)
{
	return C_f3 * (1.0 - (H_p/C_f4));
}

double f_cr(double V_TAS, double Thrust)
{
return C_fcr * n(V_TAS) * Thrust;
}

double f_ap_ld(double V_TAS, double Thrust, double H_p)
{
	return std::max(f_nom(V_TAS, Thrust), f_min(H_p));
}

//********************************************
// Maximum Climb & Take-Off Thrust [Jet Only]
//********************************************

double Thrust_max_climb(double H_p)
{
	return C_T_c1 * (1.0 - (H_p / C_T_c2) + C_T_c3 * H_p * H_p);  // (3.7-1)
}

double Thrust_max_cruise(double H_p)
{
	return C_T_CR * Thrust_max_climb(H_p);                        // (3.7-8)
}

double Descent_Thrust(double H_p, std::string Configuration)
{
	if (H_p>H_p_des){
		return C_T_des_high * Thrust_max_climb(H_p);        // (3.7-9)
	}
	else{
		if(Configuration == "Cruise")
		{
			return C_T_des_low * Thrust_max_climb(H_p); // (3.7-10)
		}
		else if(Configuration == "Approach")
		{
			return C_T_des_app * Thrust_max_climb(H_p); // (3.7-11)
		}
		else if(Configuration == "Landing")
		{
			return C_T_des_ld * Thrust_max_climb(H_p); // (3.7-12)
		}
		else
		{   	std::cout<<"\n"<<"****************************************************************"<<std::endl;
			std::cout<<"Error: Descent_Trust std::string Configuration: INVALID ARGUMENT"<<std::endl;
                        std::cout<<"Configuration can take values: Cruise, Approach or Landing"<<std::endl;
			std::cout<<"****************************************************************"<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
}

//********************************************
//           Atmosphere Model 
//********************************************

double Air_Temperature(double H_p)
{
	if (H_p<=H_p_trop)                         
	{
		return To + beta_T * H_p;          //(3.1-13)
	}
	else
	{
		return T_ISA_trop;                     //(3.1-16) 		
	}
}

double Air_Pressure(double H_p)
{
	double T = Air_Temperature(H_p);
	if (H_p<=H_p_trop)                         
	{
		return Po*std::pow((T/To), -(g/(beta_T * R)));          //(3.1-18)
	}
	else
	{
		return P_trop*exp(-(g/(R*T_ISA_trop))*(H_p-H_p_trop));  //(3.1-20)
	}
}

double Air_Density(double H_p)
{
	double T = Air_Temperature(H_p);	
	return (Air_Pressure(H_p)/(R * T));//*1.0E-06;
}

double Speed_of_Sound(double H_p)
{
	return sqrt(k * R * Air_Temperature(H_p));
}

//********************************************
// Point Mass Aircraft Model (North East UP)
//********************************************
const double S = 122.6E-06;         // [km^2]
const double C_L0 = 0.1;
const double C_La = 0.2567;
const double C_D0 = 0.038;
const double C_D2 = 0.0419;
const double Eta_max = 1.2;
const double Eta_min = 0.8;
//********************************************

double Lift(double H_p, double V_TAS, double alpha)
{
	double C_L = C_L0 + C_La*alpha;
	return (Air_Density(H_p)/2.0)*S*V_TAS*V_TAS*C_L;
}

double Drag(double H_p, double V_TAS, double alpha)
{
	double C_L = C_L0 + C_La*alpha;
	double C_D = C_D0 + C_D2*C_L*C_L;
	return  (Air_Density(H_p)/2.0)*S*V_TAS*V_TAS*C_D;
}

// Load Factor - Eta
double Eta(double H_p, double V_TAS, double alpha,double mass)
{
	return  Lift(H_p, V_TAS, alpha)/(mass*g);
}

//**************************************************************
//                  Aircraft Models
//**************************************************************

void Aircraft_Model_Climb(double *xd, double *u, double *rhs, double *p)
{
	double Tmax = 0.14104;
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);

	double Thrust = Tmax * u[1];
	double fr = f_nom(xd[2], Thrust);

	rhs[0] = p[0]*xd[2]*cos(xd[3]); 
	rhs[1] = p[0]*xd[2]*sin(xd[3]);
	rhs[2] = p[0]*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = p[0]*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -p[0]*fr;
}

void Aircraft_Model_Cruise(double *xd, double *u, double *rhs, double *p)
{
        double Tmax = 0.14104;
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);

	double Thrust = C_T_CR * Tmax * u[1];
	double fr = f_cr(xd[2], Thrust);

	rhs[0] = (p[1]-p[0])*xd[2]*cos(xd[3]); 
	rhs[1] = (p[1]-p[0])*xd[2]*sin(xd[3]);
	rhs[2] = (p[1]-p[0])*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = (p[1]-p[0])*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -(p[1]-p[0])*fr;

}

void Aircraft_Model_Descent(double *xd, double *u, double *rhs, double *p)
{
        double Tmax = 0.14104;	
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);
	
	double Thrust = Tmax * u[1];
	double fr = f_nom(xd[2], Thrust);

	rhs[0] = (p[2]-p[1])*xd[2]*cos(xd[3]); 
	rhs[1] = (p[2]-p[1])*xd[2]*sin(xd[3]);
	rhs[2] = (p[2]-p[1])*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = (p[2]-p[1])*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -(p[2]-p[1])*fr;

}


//**************************************************************
//                  Aircraft Models - Thrust
//**************************************************************


void Aircraft_Model_Climb_Thrust(double *xd, double *u, double *rhs, double *p)
{
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);
	
	double Tmax = Thrust_max_climb(xd[1]);
	double Thrust = Tmax * u[1];
	double fr = f_nom(xd[2], Thrust);


	rhs[0] = p[0]*xd[2]*cos(xd[3]); 
	rhs[1] = p[0]*xd[2]*sin(xd[3]);
	rhs[2] = p[0]*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = p[0]*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -p[0]*fr;

}

void Aircraft_Model_Cruise_Thrust(double *xd, double *u, double *rhs, double *p)
{
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);
	
	double Tmax = Thrust_max_cruise(xd[1]);
	double Thrust = Tmax * u[1];
	double fr = f_cr(xd[2], Thrust);

	rhs[0] = (p[1]-p[0])*xd[2]*cos(xd[3]); 
	rhs[1] = (p[1]-p[0])*xd[2]*sin(xd[3]);
	rhs[2] = (p[1]-p[0])*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = (p[1]-p[0])*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -(p[1]-p[0])*fr;



}


void Aircraft_Model_Descent_Thrust_high(double *xd, double *u, double *rhs, double *p)
{
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);
	
	double Tmax = Descent_Thrust(xd[1],"Cruise");

	double Thrust = Tmax * u[1];
	double fr = f_nom(xd[2], Thrust);

	rhs[0] = (p[2]-p[1])*xd[2]*cos(xd[3]); 
	rhs[1] = (p[2]-p[1])*xd[2]*sin(xd[3]);
	rhs[2] = (p[2]-p[1])*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = (p[2]-p[1])*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -(p[2]-p[1])*fr;

}

void Aircraft_Model_Descent_Thrust_low(double *xd, double *u, double *rhs, double *p)
{
	double L = Lift(xd[1],xd[2],u[0]);
	double D = Drag(xd[1],xd[2],u[0]);
	
	double Tmax = Descent_Thrust(xd[1],"Cruise");

	double Thrust = Tmax * u[1];
	double fr = f_nom(xd[2], Thrust);

	rhs[0] = (p[3]-p[2])*xd[2]*cos(xd[3]); 
	rhs[1] = (p[3]-p[2])*xd[2]*sin(xd[3]);
	rhs[2] = (p[3]-p[2])*(1.0/xd[4])*(Thrust - D)-g*sin(xd[3]);
	rhs[3] = (p[3]-p[2])*(1.0/(xd[4]*xd[2]))*(L-g*xd[4]*cos(xd[3]));
	rhs[4] = -(p[3]-p[2])*fr;

}

#endif
