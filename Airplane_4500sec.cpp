#include <cmath>
#include <iostream>
#include <iomanip>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include "def_usrmod.hpp"
#include <fstream>

#include "Models_v6.cpp"    

#define  NMOS   3  /* Number of phases (Model Stages) */
#define  NP     3  /* Number of parameters */
#define  NRC    0  /* Number of coupled constraints */
#define  NRCE   0  /* Number of coupled equality constraints */

#define  NXD    5  /* Number of differential states */
#define  NXA    0  /* Number of algebraic states */
#define  NU     2  /* Number of controls */
#define  NPR    0  /* Number of local parameters */

//Climb
#define  NRD_CL_S  7  /* Number of constraints at the start point */
#define  NRDE_CL_S 4  /* Number of equality constraints at the start points */

#define  NRD_CL_E  3  /* Number of constraints at the end point */
#define  NRDE_CL_E 1  /* Number of equality constraints at the end point */

//Cruise
#define  NRD_CR_E  3  /* Number of constraints at the end point */
#define  NRDE_CR_E 1  /* Number of equality constraints at the end point */

//Descent
#define  NRD_DE_E  6  /* Number of constraints at the end point */
#define  NRDE_DE_E 3  /* Number of equality constraints at the end point */

//Interior Points
#define NRD_I   4  /* Number of constraints at interior points */
#define NRDE_I  0  /* Number of equality constraints at interior points */

// Name of .csv file to save Optimisation results 

std::string result_file="RES/Result_4500sec.csv";
std::string param_file="RES/Parameters_4500sec.csv";
// Objective Function Weights 

const double q_CL_d = 0.0001;

const double q_CR_d = 0.0001;

const double q_DE_d = 0.0001;
const double q_DE_m = 5.0;

const double q_DE_t = 0.0;

double M_TO =70.0;

double V_TO = 0.074594444;//268.54/3600.0;0,074594444
double V_LD = 0.069833333;//251.4/3600.0;0,069833333

double X_LD   = 1000.0;
double H_p_LD = 0.0;

double X_TO   = 0.0;
double H_p_TO = 0.0;

//************************************************************************
//           Lagrange Objective Function Fuel Consumption - Climb
//************************************************************************
static void lfcn_climb_fuel(double *t, double *xd, double *xa, double *u,
  double *p, double *lval, double *rwh, long *iwh, InfoPtr *info) {

	*lval = q_CL_d*p[0]*u[1]*u[1]; 
}

//************************************************************************
//           Lagrange Objective Function Fuel Consumption - Cruise
//************************************************************************
static void lfcn_cruise_fuel(double *t, double *xd, double *xa, double *u,
  double *p, double *lval, double *rwh, long *iwh, InfoPtr *info) {

	*lval = q_CR_d*(p[1]-p[0])*u[1]*u[1];

}

//************************************************************************
//           Lagrange Objective Function Fuel Consumption - Descent
//************************************************************************
static void lfcn_descent_fuel(double *t, double *xd, double *xa, double *u,
  double *p, double *lval, double *rwh, long *iwh, InfoPtr *info) {

	*lval = q_DE_d*(p[2]-p[1])*u[1]*u[1];

}
//***********************************************************
//           Mayer Objective Function End Time - Descent
//***********************************************************
static void mfcn_descent_end_time(double *ts, double *xd, double *xa, double *p, double *pr,
    double *mval,  long *dpnd, InfoPtr *info) {
  if (*dpnd) {
      *dpnd = MFCN_DPND(0, *xd, 0, 0, 0);   //*ts, *sd, *sa, *u, *p, *pr
      return;
  }

  *mval = - q_DE_m * xd[4];
}

//***************************************************
//alpha --> u[0], delta --> u[1] ,mu --> u[2]
//***************************************************
//            Right Hand Side of Climb Stage
//***************************************************
static void ffcn_climb(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
	Aircraft_Model_Climb(xd, u, rhs, p);
}
//***************************************************
//            Right Hand Side of Cruise Stage
//***************************************************

static void ffcn_cruise(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
	Aircraft_Model_Cruise(xd, u, rhs, p);
}
//***************************************************
//            Right Hand Side of Descent Stage
//***************************************************
static void ffcn_descent(double *t, double *xd, double *xa, double *u,
  double *p, double *rhs, double *rwh, long *iwh, InfoPtr *info)
{
	Aircraft_Model_Descent(xd, u, rhs, p);
}
//***************************************************
//            Constraints at Start point - Climb
//***************************************************
static void rdfcn_climb_s(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
    return;
  }
  res[0] = sd[0]-X_TO;
  res[1] = sd[1]-H_p_TO;
  res[2] = (sd[2]-V_TO);
  res[3] = (sd[4]-M_TO);
  res[4] = sd[3];
  res[5] = (Eta(sd[1], sd[2], u[0], sd[4])-Eta_min);
  res[6] = (Eta_max-Eta(sd[1], sd[2], u[0], sd[4]));
}
//***************************************************
//            Interior Point Constraints 
//***************************************************
static void rdfcn_i(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
	if (*dpnd) { *dpnd = RFCN_DPND(0, *sd, 0, *u, *p, 0); return; }

  res[0] = (Eta(sd[1], sd[2], u[0], sd[4])-Eta_min);
  res[1] = (Eta_max-Eta(sd[1], sd[2], u[0], sd[4]));
  res[2] = p[1]-p[0];
  res[3] = p[2]-p[1];

}
//***************************************************
//            Constraints at End point - Climb
//***************************************************
static void rdfcn_climb_e(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
    return;
  }
  //res[0] = sd[1] - H_p_CR;
  res[0] = sd[3];
  res[1] = (Eta(sd[1], sd[2], u[0], sd[4])-Eta_min);
  res[2] = (Eta_max-Eta(sd[1], sd[2], u[0], sd[4]));
}
//***************************************************
//            Constraints at End point - Cruise
//***************************************************
static void rdfcn_cruise_e(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
    return;
  }
  //res[0] = sd[1] - H_p_CR;
  res[0] = sd[3];
  res[1] = (Eta(sd[1], sd[2], u[0], sd[4])-Eta_min);
  res[2] = (Eta_max-Eta(sd[1], sd[2], u[0], sd[4]));
}
//***************************************************
//            Constraints at End point - Descent
//***************************************************
static void rdfcn_descent_e(double *ts, double *sd, double *sa, double *u,
  double *p, double *pr, double *res, long *dpnd, InfoPtr *info)
{
  if (*dpnd) {
    *dpnd = RFCN_DPND(0, *sd, 0, *u, 0, 0);
    return;
  }
  res[0] = sd[0]-X_LD;
  res[1] = sd[1]-H_p_LD;
  res[2] = sd[2]-V_LD;
  res[3] = -sd[3];
  res[4] = (Eta(sd[1], sd[2], u[0], sd[4])-Eta_min);
  res[5] = (Eta_max-Eta(sd[1], sd[2], u[0], sd[4]));
}
//***************************************************
//            Save Results
//***************************************************
std::vector<double> t_values;
std::vector<double> eta_values;
std::vector<std::vector<double>> sd_values;
std::vector<std::vector<double>> u_values;

static void data_out( double *t, double *sd, double *sa, double *u, double *p, double *rwh, long *iwh, InfoPtr *info ) {
	if (*t == 0.) {
		std::ofstream csv_stream;
		std::ofstream param_stream;

		csv_stream.open(result_file.c_str(), std::ios_base::trunc);
		param_stream.open(param_file.c_str(), std::ios_base::trunc);

		if (t_values.size() > 0) {
			param_stream << p[0] << "\t"<< p[1] << "\t"<< p[2] << std::endl;
			for (unsigned int i = 0; i < t_values.size(); i ++) {
				csv_stream << t_values[i] << "\t";
				for (unsigned int j = 0; j < sd_values[i].size(); j++) {
					csv_stream << sd_values[i][j];
					if (j < sd_values[i].size() -1 )
						csv_stream << "\t";
				}
                                csv_stream << "\t";
                                for (unsigned int j = 0; j < u_values[i].size(); j++) {
					csv_stream << u_values[i][j];
					if (j < u_values[i].size() -1 )
						csv_stream << "\t";
				}
//				csv_stream << std::endl;
				csv_stream << "\t" <<eta_values[i] << std::endl;
			}
		}

		t_values.clear();
		sd_values.clear();
		u_values.clear();
                eta_values.clear();

		csv_stream.close();
		param_stream.close();
	}

	t_values.push_back(t[0]);

        eta_values.push_back(Eta(sd[1], sd[2], u[0], sd[4]));

	std::vector<double> sd_vec(NXD);
	for (unsigned i = 0; i < NXD; i++)
		sd_vec[i] = sd[i];
	sd_values.push_back(sd_vec);

	std::vector<double> u_vec(NU);
	for (unsigned i = 0; i < NU; i++)
		u_vec[i] = u[i];
	u_values.push_back(u_vec);
}

void mout
(
  long   *imos,      ///< index of model stage (I)
  long   *imsn,      ///< index of m.s. node on current model stage (I)
  double *ts,        ///< time at m.s. node (I)
  double *te,        ///< time at end of m.s. interval (I)
  double *sd,        ///< differential states at m.s. node (I)
  double *sa,        ///< algebraic states at m.s. node (I)
  double *u,         ///< controls at m.s. node (I)
  double *udot,      ///< control slopes at m.s. node (I)
  double *ue,        ///< controls at end of m.s. interval (I)
  double *uedot,     ///< control slopes at end of m.s. interval (I)
  double *p,         ///< global model parameters (I)
  double *pr,        ///< local i.p.c. parameters (I)
  double *ccxd,
  double *mul_ccxd,  ///< multipliers of continuity conditions (I)
#if defined(PRSQP) || defined(EXTPRSQP)
  double *ares,
  double *mul_ares,
#endif
  double *rd,
  double *mul_rd,    ///< multipliers of decoupled i.p.c. (I)
  double *rc,
  double *mul_rc,    ///< multipliers of coupled i.p.c. (I)
  double *obj,
  double *rwh,       ///< real work array (I)
  long   *iwh        ///< integer work array (I)
)
{
    InfoPtr info(0, *imos, *imsn);
    data_out( ts, sd, sa, u, p, rwh, iwh, &info);
}
//***************************************************
//            Model Definition
//***************************************************
extern "C" void def_model(void);
void def_model(void)
{
	/* Define problem dimensions */
	def_mdims(NMOS, NP, NRC, NRCE);

// Define the Climb Phase 
	
		def_mstage(
			0,
			NXD, NXA, NU,
			NULL, lfcn_climb_fuel, // mayer term, lagrange term
                        0, 0, 0, NULL, ffcn_climb, NULL,
			NULL, NULL
			);

// Define the Cruise Phase 
	
		def_mstage(
			1,
			NXD, NXA, NU,
			NULL, lfcn_cruise_fuel, // mayer term, lagrange term
                        0, 0, 0, NULL, ffcn_cruise, NULL,
			NULL, NULL
			);

// Define the Descent Phase 
	
		def_mstage(
			2,
			NXD, NXA, NU,
			mfcn_descent_end_time, lfcn_descent_fuel, // mayer term, lagrange term
                        0, 0, 0, NULL, ffcn_descent, NULL,
			NULL, NULL
			);

//********************************************************************
//        Define constraints at the Start, Interior, End point 
//********************************************************************
	def_mpc(0, "Start Point", NPR, NRD_CL_S, NRDE_CL_S, rdfcn_climb_s, NULL);
        def_mpc(0, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);

	def_mpc(1, "Start Point", NPR, NRD_CL_E, NRDE_CL_E, rdfcn_climb_e, NULL);
        def_mpc(1, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);

	def_mpc(2, "Start Point", NPR, NRD_CR_E, NRDE_CR_E, rdfcn_cruise_e, NULL);
        def_mpc(2, "Interior Point", NPR, NRD_I, NRDE_I, rdfcn_i, NULL);
	def_mpc(2, "End Point", NPR, NRD_DE_E, NRDE_DE_E, rdfcn_descent_e, NULL);

        def_mio(NULL ,mout,data_out);
}

