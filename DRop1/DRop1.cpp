
#include <stdio.h>
#include <math.h>
#include<cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif
/* Problem Constants */
#define H2O  0
#define N2 1
#define INTERFACE 2
#define NEQ   551
#define TOUT  0.2
#define WIDTH 10
#define N RCONST(NEQ)
#define N_d 50
#define l0 RCONST(0.01)
#define L RCONST(5.5)

#define lambda_H2O RCONST(0.56e-3)
#define lambda_N2 RCONST(0.05e-3)
#define Cp_H20 RCONST(4180.0)
#define Cp_N2 RCONST(2000)
#define rho_H20 RCONST(0.98e-6)
#define rho_N2 RCONST(4.0e-9)
#define T0_H2O RCONST(300)
#define T0_N2 RCONST(1500)
#define T_vep RCONST(373)
#define Ld RCONST(2258.2e3)
#define R RCONST(8.31)
#define p0 1.2
#define V0 0
#define tau RCONST(1.0e-4)
#define tau0 RCONST(1.0e-4)
#define P0 RCONST(1.e5)
#define M_N2 RCONST(0.018)

# define M_PI 3.14159265358979323846
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
/* Macro to define dense matrix elements, indexed from 1. */

#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)

typedef struct {
    realtype x[N - 1];
    realtype p;
    double p_new;
    int inter;
    realtype t;
    realtype tout;
    bool step;
    void* mem;
} *UserData;


/* Prototypes of functions called by IDA */

int resHEAT(realtype tres, N_Vector yy, N_Vector yp,
    N_Vector resval, void* user_data);
int resHEAT_inital(realtype tres, N_Vector yy, N_Vector yp,
    N_Vector resval, void* user_data);
int resEVEPORATION(realtype tres, N_Vector yy, N_Vector yp,
    N_Vector resval, void* user_data);


/*int jacrob(realtype tt, realtype cj,
    N_Vector yy, N_Vector yp, N_Vector resvec,
    SUNMatrix JJ, void* user_data,
    N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);*/

    /* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);
static void PrintOutput(void* mem, realtype t, N_Vector y);
static void PrintRootInfo(int root_f1, int root_f2);
static int check_retval(void* returnvalue, const char* funcname, int opt);
static int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol);
/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

realtype lambda(int phase) {
    if (phase == 0) {
        return lambda_H2O;
    }
    if (phase == 1) {
        return lambda_N2;
    }
    if (phase == 2) {
        return 0.5 * (lambda_H2O + lambda_N2);
    }
}

realtype rho(int phase, double T) {
    if (phase == 0) {
        return rho_H20;
    }
    else {
        return(P0 * M_N2 / (R * T)) * 1.e-9;;
    }
}
realtype Cp(int phase) {
    if (phase == 0) {
        return Cp_H20;
    }
    else {
        return Cp_N2;
    }
}

double Func_therm(double riM, double ri, double riP, double Tvall, double Tval, double Tvalr, int phaseL, int phaseR) {
    return 1 / pow(ri, 2) * (2 / (riP - riM)) *
        (lambda(phaseR) * pow(0.5 * (riP + ri), 2) * (Tvalr - Tval) / (riP - ri) -
            lambda(phaseL) * pow(0.5 * (ri + riM), 2) * (Tval - Tvall) / (ri - riM));
}
double Func_therm_inter_P(double riM, double ri, double riP, double Tvall, double Tval, double Tvalr, double TvalrP, int phaseL, int phaseR, double h, double p) {
    return 1 / pow(ri, 2) * (1 / (0.5 * (riP + ri) - riM)) *
        (lambda(phaseL) * pow(0.5 * (riP + ri), 2) * (Tvalr - Tval) / (riP - ri) -
            pow(riM, 2) * lambda(phaseR) * (((2 * p - 7) / ((p - 3) * (p - 4) * h)) * Tvall + ((p - 4) / ((p - 3) * h)) * Tvalr + ((p - 3) / (h * (4 - p))) * TvalrP));
}
double Func_therm_inter_M(double riM, double ri, double riP, double TvallM, double Tvall, double Tval, double Tvalr, int phaseL, int phaseR, double h, double p) {
    return  1 / pow(ri, 2) * (1 / (riP - 0.5 * (riM + ri))) *
        (pow(riP, 2) * lambda(phaseR) * (p / ((p + 1) * h) * TvallM - (p + 1) / (p * h) * Tvall + Tvalr * (2 * p + 1) / ((p + 1) * p * h)) -
            lambda(phaseL) * pow(0.5 * (ri + riM), 2) * (Tval - Tvall) / (ri - riM));
}
double Func_therm_inter(double TvallM, double Tvall, double Tval, double Tvalr, double TvalrP, int phaseL, int phaseR, double h, double p) {
    return lambda(phaseL) * (p / (p + 1) * TvallM - (p + 1) / p * Tvall + Tval * (2 * p + 1) / ((p + 1) * p)) / h -
        lambda(phaseR) * (((2 * p - 7) / ((p - 3) * (p - 4))) * Tval + ((p - 4) / (p - 3)) * Tvalr + ((p - 3) / (4 - p)) * TvalrP) / h;
}


double Func_therm_inter_evoparation_P(double riM, double ri, double riP, double Tvall, double Tval, double Tvalr, double TvalrP, double dot_M, int phaseL, int phaseR, double h, double p) {
    return 1 / pow(ri, 2) * (1 / (0.5 * (riP + ri) - riM)) *
        (lambda(phaseL) * pow(0.5 * (riP + ri), 2) * (Tvalr - Tval) / (riP - ri) -
            pow(riM, 2) * lambda(phaseR) * (((2 * p - 7) / ((p - 3) * (p - 4) * h)) * Tvall + ((p - 4) / ((p - 3) * h)) * Tvalr + ((p - 3) / (h * (4 - p))) * TvalrP)) - dot_M * Cp(phaseR) * (Tval - Tvall) / ((ri - riM) * pow(ri, 2));
}
double Func_therm_evoparation_inter(double TvallM, double Tvall, double Tval, double Tvalr, double TvalrP, double dot_M, double ri, int phaseL, int phaseR, double h, double p) {
    return lambda(phaseL) * (p / (p + 1) * TvallM - (p + 1) / p * Tvall + Tval * (2 * p + 1) / ((p + 1) * p)) / h -
        lambda(phaseR) * (((2 * p - 7) / ((p - 3) * (p - 4))) * Tval + ((p - 4) / (p - 3)) * Tvalr + ((p - 3) / (4 - p)) * TvalrP) / h + Ld * dot_M / (ri * ri);
}
double Func_therm_evoparation(double riM, double ri, double riP, double Tvall, double Tval, double Tvalr, double dot_M, int phaseL, int phaseR) {
    return 1 / pow(ri, 2) * (2 / (riP - riM)) *
        (lambda(phaseR) * pow(0.5 * (riP + ri), 2) * (Tvalr - Tval) / (riP - ri) -
            lambda(phaseL) * pow(0.5 * (ri + riM), 2) * (Tval - Tvall) / (ri - riM)) - dot_M * Cp(phaseR) * (Tval - Tvall) / ((ri - riM) * pow(ri, 2));
}

double get_q(double Tval, double Tvalr, double TvalrP, int phaseR, double h, double p) {
    return lambda(phaseR) * (((2 * p - 7) / ((p - 3) * (p - 4))) * Tval + ((p - 4) / (p - 3)) * Tvalr + ((p - 3) / (4 - p)) * TvalrP) / h;
}
double get_p(double t, double dot_M, double ri, double h, void* user_data, void* mem) {
    UserData data;
    data = (UserData)user_data;
    double p;
    int retval;
    double V;
    double dt;
    double tout;
    double step;
    double tprev;
    int inter;
    inter = data->inter;
    tprev = data->t;
    retval = IDAGetCurrentStep(mem, &step);
    check_retval(&retval, "IDAGetCurrentStep", 1);
    tout = data->tout;
    dt = t - tprev;
    p = data->p;
    data->t = t;
    if (abs(dt) > 1.0e-17) {
        data->p = data->p_new;
        data->t = t;
        p = data->p;
    }
    p = data->p;
    V = dot_M * 1 / (rho_H20 * ri * ri * h);
    p -= step * V;
    if (abs(p - 1.0) < 0.001) {
        p = 1.999;
        inter = inter - 1;
        data->inter = inter;
    }
    data->p_new = p;
    data->tout = tout;
    return p;

}
std::vector<realtype> mesh() {
    int inter = N_d - 1;
    int cout = 0;
    realtype h = L / (N - 1);
    std::vector<realtype> X;
    for (int i = 0; i < N - 1; i++) {
        X.push_back(i * h + 0.5 * h);
        //std::cout << R.at(i) << "   " << i << std::endl;
    }
    return X;
}
void out_T(double tout, double* Tval, void* user_data) {
    std::string name;
    double ri;
    double p;
    double h;
    UserData data;
    data = (UserData)user_data;
    int inter;
    inter = data->inter;
    p = data->p;
    h = data->x[inter - 3] - data->x[inter - 4];
    name = "output" + std::to_string(tout) + std::to_string(N - 1) + ".txt";
    std::ofstream outputFile;
    outputFile.open(name);
    for (int i = 0; i < NEQ; i++) {
        if (i < inter + 1) {
            ri = data->x[i];
            outputFile << ri << "," << Tval[i] << std::endl;
        }
        /*else if (i == inter + 1) {
            ri = data->x[inter - 1] + (p)*h;
            outputFile << ri << "," << Tval[i] << std::endl;
            i++;
        }*/
        if (i == inter + 1) {
            ri = data->x[inter - 1] + (p)*h;
            outputFile << ri << "," << T_vep << std::endl;
        }
        if (i > inter + 1) {

            ri = data->x[i - 1];
            outputFile << ri << "," << Tval[i] << std::endl;
        }

        /* ri = data->x[i];
         outputFile << ri << "," << Tval[i + 1] << std::endl;
         std::cout<< ri << "," << Tval[i + 1] << std::endl;*/
    }
    outputFile.close();
}
void output(double tout, double* T, void* user_data) {
    std::string name;
    double h, ri, p;
    int inter;
    UserData data;
    data = (UserData)user_data;
    p = data->p;
    inter = data->inter;
    h = data->x[inter - 1] - data->x[inter - 2];
    name = "out\\output" + std::to_string(tout) + ".dat";
    std::ofstream outputFile;
    outputFile.open(name);
    outputFile << "TITLE=\"" << "Graphics" << "\"" << std::endl;
    outputFile << R"(VARIABLES= "r", "T")" << std::endl;
    for (int i = 0; i < NEQ; i++) {
        if (i < inter + 1) {
            ri = data->x[i];
            outputFile << ri << " " << T[i] << std::endl;

        }

        if (i == inter + 1) {
            ri = data->x[inter - 1] + (p)*h;
            outputFile << ri << " " << T[inter + 1] << std::endl;
        }
        if (i > inter + 1) {
            ri = data->x[i - 1];
            outputFile << ri << " " << T[i] << std::endl;
        }
        outputFile.close();

    }
}
void out_T1(double tout, double* Tval, void* user_data) {
    std::string name;
    double ri;
    double p;
    double h;
    UserData data;
    data = (UserData)user_data;
    int inter;
    inter = data->inter;
    p = data->p;
    h = data->x[inter - 3] - data->x[inter - 4];
    name = "output" + std::to_string(tout) + std::to_string(N - 1) + ".txt";
    std::ofstream outputFile;
    outputFile.open(name);
    for (int i = 0; i < NEQ - 1; i++) {
        ri = data->x[i];
        outputFile << ri << "," << Tval[i] << std::endl;

        /* ri = data->x[i];
         outputFile << ri << "," << Tval[i + 1] << std::endl;
         std::cout<< ri << "," << Tval[i + 1] << std::endl;*/
    }
    outputFile.close();
}
void out_Inter(double tout, double* Tval, void* user_data) {
    std::string name;
    double r_inter, h, p;
    int inter;
    UserData data;
    data = (UserData)user_data;
    inter = data->inter;
    name = "Inter" + std::to_string(tout) + ".txt";
    h = data->x[inter - 1] - data->x[inter - 2];
    p = data->p;
    r_inter = data->x[inter - 1] + p * h;
    std::cout << r_inter << "***" << std::endl;
    std::ofstream outputFile;
    outputFile.open(name);
    outputFile << r_inter << "," << Tval[inter + 1] << std::endl;
    outputFile.close();
}
//double get_v(double dot_M, double ri,double h) {
//    return dot_M /(4*M_PI* rho(H2O) * ri * dt * h)
//}
int inital_T(double(&Init_T)[NEQ - 1], void* user_data) {
    UserData data;
    data = (UserData)user_data;
    void* mem;
    N_Vector T, Tp, avtol;
    realtype rtol, * Tval, * Tpval, * atval;
    realtype t0, tout1, tout, tret, h, ri, riM, riP, v;
    int iout, retval, retvalr;
    int inter = N_d - 1;
    int rootsfound[2];
    SUNMatrix B;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    h = L / (N - 1);
    mem = NULL;
    T = Tp = avtol = NULL;
    Tval = Tpval = atval = NULL;
    B = NULL;
    LS = NULL;
    NLS = NULL;
    // std::cout << "WIDTH + N_d - 2 = " << R_vec.at(WIDTH + N_d - 2) << "\n";
    /*data->inter = inter;*/
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);
    T = N_VNew_Serial(NEQ - 1, ctx);
    if (check_retval((void*)T, "N_VNew_Serial", 0)) return(1);
    Tp = N_VClone(T);
    if (check_retval((void*)Tp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(T);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    Tval = N_VGetArrayPointer(T);
    // std::cout << "Ok" << "\n";
     //alpha = 300.0 / (data->x[WIDTH + N_d - 2 + l] - data->x[WIDTH + N_d - 2 - l]);
    //Тут размер массива x и Т разных размерностей
    for (int i = 0; i < NEQ - 1; i++) {
        if (i < inter + 1) {
            Tval[i] = T0_H2O;
        }
        else {
            Tval[i] = T0_N2;
        }
    }
    h = data->x[inter] - data->x[inter - 1];
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);
    for (int i = 0; i < NEQ - 1; i++) {
        atval[i] = RCONST(1.0e-3);
    }
    Tpval = N_VGetArrayPointer(Tp);
    ri = data->x[0];
    riM = 0.0;
    riP = data->x[1];
    Tpval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[0]);
    for (int i = 1; i < inter; i++) {
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[i]);// - функция справа от производных 
        std::cout << "F = " << riP << " i=  " << i << std::endl;
    }
    Tpval[inter] = Func_therm(riM, ri, riP, Tval[inter - 1], Tval[inter], Tval[inter + 1], H2O, N2) / Cp(H2O) / rho(H2O, Tval[inter]);
    Tpval[inter + 1] = Func_therm(riM, ri, riP, Tval[inter], Tval[inter + 1], Tval[inter + 2], H2O, N2) / Cp(N2) / rho(N2, Tval[inter]);
    for (int i = inter + 2; i < NEQ - 2; i++) {
        ri = data->x[i - 1];
        riP = data->x[i - 1] + h;
        riM = data->x[i - 2];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], N2, N2) / Cp(N2) / rho(N2, Tval[i]);// - функция справа от производных
    }
    ri = data->x[NEQ - 2];
    riP = data->x[NEQ - 2] + h;
    riM = data->x[NEQ - 3];
    Tpval[NEQ - 2] = Func_therm(riM, ri, riP, Tval[NEQ - 3], Tval[NEQ - 2], Tval[NEQ - 2], N2, N2) / Cp(N2) / rho(N2, Tval[NEQ - 2]);
    /* Integration limits */
    t0 = ZERO;
    tout1 = tau0;

    PrintHeader(rtol, avtol, T);

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);
    retval = IDAInit(mem, resHEAT_inital, t0, T, Tp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    data->mem = mem;

    /* Create dense SUNMatrix for use in linear solves */
    B = SUNDenseMatrix(NEQ - 1, NEQ - 1, ctx);
    if (check_retval((void*)B, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(T, B, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, B);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine */
   // retval = IDASetJacFn(mem, jacrob);
    //if (check_retval(&retval, "IDASetJacFn", 1)) return(1);

    /* Create Newton SUNNonlinearSolver object. IDA uses a
     * Newton SUNNonlinearSolver by default, so it is unecessary
     * to create it and attach it. It is done in this example code
     * solely for demonstration purposes. */
    NLS = SUNNonlinSol_Newton(T, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */
    ri = data->x[inter];
    riP = data->x[inter + 1];
    riM = data->x[inter - 1];
    iout = 0; tout = tout1;
    data->tout = 0;
    while (1) {
        //std::cout << "inter  " << data->x[inter - 1] + p * h << std::endl;
        retval = IDASolve(mem, tout, &tret, T, Tp, IDA_NORMAL);

        if (check_retval(&retval, "IDASolve", 1)) return(1);

        if (retval == IDA_ROOT_RETURN) {
            retvalr = IDAGetRootInfo(mem, rootsfound);
            check_retval(&retvalr, "IDAGetRootInfo", 1);
            PrintRootInfo(rootsfound[0], rootsfound[1]);
        }

        if (retval == IDA_SUCCESS) {
            iout++;
            data->t = tout;
            tout += tau0;
            iout++;
            //output(tout, Tval, data);
        }
        // output(tout, Tval, data);
        if (iout == 100000) break;
    }
    for (int i = 0; i < NEQ - 1; i++) {
        Init_T[i] = Tval[i];
    }

    std::cout << "t =   " << tout << std::endl;
    printf("\n\n");
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(B);
    N_VDestroy(avtol);
    N_VDestroy(T);
    N_VDestroy(Tp);
    SUNContext_Free(&ctx);


}
int heat_inter(double(&inter_T)[NEQ], void* user_data) {
    UserData data;
    data = (UserData)user_data;
    void* mem;
    N_Vector T, Tp, avtol;
    realtype rtol, * Tval, * Tpval, * atval;
    realtype t0, tout1, tret, tout, h, ri, riM, riP, p;
    int iout, retval, retvalr;
    int inter = N_d - 1;
    int rootsfound[2];
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    h = L / (N - 1);
    mem = NULL;
    T = Tp = avtol = NULL;
    Tval = Tpval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
    p = data->p;
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);
    T = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)T, "N_VNew_Serial", 0)) return(1);
    Tp = N_VClone(T);
    if (check_retval((void*)Tp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(T);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    Tval = N_VGetArrayPointer(T);
    for (int i = 0; i < NEQ; i++) {
        /*Tval[i] = inital_T(ri);*/
        if (i < inter + 1) {
            Tval[i] = T0_H2O;
        }
        if (i > inter + 1) {
            Tval[i] = T0_N2;
        }
    }

    Tval[inter + 1] = 370.61;
    //(lambda(H2O) * T0_H2O + lambda(N2) * Tval[inter + 2]) / (lambda(N2) + lambda(H2O));
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);
    for (int i = 0; i < NEQ; i++) {
        atval[i] = RCONST(1.0e-3);
    }
    Tpval = N_VGetArrayPointer(Tp);
    ri = data->x[0];
    riM = 0.0;
    riP = data->x[1];
    Tpval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[0]);
    for (int i = 1; i < inter; i++) {
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[i]);// - функция справа от производных 
        // std::cout << Tpval[i] << " i=  " << i << std::endl;
    }
    h = data->x[inter] - data->x[inter - 1];
    riP = data->x[inter - 1] + p * h;
    riM = data->x[inter - 1];
    ri = data->x[inter];
    Tpval[inter] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 3], Tval[inter - 2], Tval[inter - 1], Tval[inter + 1], N2, N2, h, p) / Cp(H2O) / rho(H2O, Tval[inter]);
    riM = data->x[inter - 1];
    Tpval[inter + 1] = 0;
    ri = data->x[inter + 1];
    riP = data->x[inter + 2];
    riM = data->x[inter - 1] + p * h;
    Tpval[inter + 2] = Func_therm_inter_P(riM, ri, riP, Tval[inter + 1], Tval[inter + 2], Tval[inter + 3], Tval[inter + 4], N2, N2, h, p) / Cp(N2) / rho(N2, Tval[inter + 1]);
    for (int i = inter + 2; i < NEQ - 1; i++) {
        ri = data->x[i - 1];
        riP = data->x[i - 1] + h;
        riM = data->x[i - 2];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], N2, N2) / Cp(N2) / rho(N2, Tval[i]);// - функция справа от производных
    }
    ri = data->x[NEQ - 2];
    riP = data->x[NEQ - 2] + h;
    riM = data->x[NEQ - 3];
    Tpval[NEQ - 1] = Func_therm(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], N2, N2) / Cp(N2) / rho(N2, Tval[NEQ - 1]);
    /* Integration limits */
    t0 = ZERO;
    tout1 = tau;
    PrintHeader(rtol, avtol, T);

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);
    retval = IDAInit(mem, resHEAT, t0, T, Tp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    IDASetMaxStep(mem, tau);
    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    data->mem = mem;

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(T, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine */
   // retval = IDASetJacFn(mem, jacrob);
    //if (check_retval(&retval, "IDASetJacFn", 1)) return(1);

    /* Create Newton SUNNonlinearSolver object. IDA uses a
     * Newton SUNNonlinearSolver by default, so it is unecessary
     * to create it and attach it. It is done in this example code
     * solely for demonstration purposes. */
    NLS = SUNNonlinSol_Newton(T, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */
    ri = data->x[inter];
    riP = data->x[inter + 1];
    riM = data->x[inter - 1];
    iout = 0; tout = tout1;
    data->tout = 0;
    while (Tval[inter + 1] < T_vep) {
        //std::cout << "inter  " << data->x[inter - 1] + p * h << std::endl;
        inter = data->inter;
        ri = data->x[inter - 1] + p * h;
        retval = IDASolve(mem, tout, &tret, T, Tp, IDA_NORMAL);

        if (check_retval(&retval, "IDASolve", 1)) return(1);

        if (retval == IDA_ROOT_RETURN) {
            retvalr = IDAGetRootInfo(mem, rootsfound);
            check_retval(&retvalr, "IDAGetRootInfo", 1);
            PrintRootInfo(rootsfound[0], rootsfound[1]);
        }
        if (iout % 1000 == 0) {
            out_T(tout, Tval, data);
        }
        if (retval == IDA_SUCCESS) {
            iout++;
            data->t = tout;
            tout += tau;
            //inter = data->inter;

        }
        double delta = TOUT - tout;
    }
    for (int i = 0; i < NEQ; i++) {
        inter_T[i] = Tval[i];
        //std::cout << "T "<<i <<" " << Tval[i] << std::endl;
    }
    //out_Inter(tout, Tval, data);
    std::cout << "Ti " << Tval[NEQ - 1] << std::endl;
    out_T(tout, Tval, data);
    printf("\n\n");
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(T);
    N_VDestroy(Tp);
    SUNContext_Free(&ctx);

    return(retval);


}
int evaporation(double(&inter_T)[NEQ], void* user_data) {
    UserData data;
    data = (UserData)user_data;
    void* mem;
    N_Vector T, Tp, avtol;
    realtype rtol, * Tval, * Tpval, * atval;
    realtype t0, tout1, tret, tout, h, ri, riM, riP, p;
    int iout, retval, retvalr;
    int inter = N_d - 1;
    int rootsfound[2];
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    h = L / (N - 1);
    mem = NULL;
    T = Tp = avtol = NULL;
    Tval = Tpval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
    p = data->p;
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);
    T = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)T, "N_VNew_Serial", 0)) return(1);
    Tp = N_VClone(T);
    if (check_retval((void*)Tp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(T);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    Tval = N_VGetArrayPointer(T);
    for (int i = 0; i < NEQ; i++) {
        Tval[i] = inter_T[i];
    }
    Tval[inter + 1] = 0.0;
    rtol = RCONST(1.0e-2);
    atval = N_VGetArrayPointer(avtol);
    for (int i = 0; i < NEQ; i++) {
        atval[i] = RCONST(1.0e-2);
    }
    Tpval = N_VGetArrayPointer(Tp);
    ri = data->x[0];
    riM = 0.0;
    riP = data->x[1];
    Tpval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[0]);
    for (int i = 1; i < inter; i++) {
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) / Cp(H2O) / rho(H2O, Tval[i]);// - функция справа от производных 
        // std::cout << Tpval[i] << " i=  " << i << std::endl;
    }
    h = data->x[inter] - data->x[inter - 1];
    riP = data->x[inter - 1] + p * h;
    riM = data->x[inter - 1];
    ri = data->x[inter];
    Tpval[inter] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 3], Tval[inter - 2], Tval[inter - 1], Tval[inter + 1], N2, N2, h, p) / Cp(H2O) / rho(H2O, Tval[inter]);
    riM = data->x[inter - 1];
    Tpval[inter + 1] = 0;// - функция справа от производных 
    ri = data->x[inter + 1];
    riP = data->x[inter + 2];
    riM = data->x[inter - 1] + p * h;
    Tpval[inter + 2] = Func_therm_inter_P(riM, ri, riP, Tval[inter + 1], Tval[inter + 2], Tval[inter + 3], Tval[inter + 4], N2, N2, h, p) / Cp(N2) / rho(N2, Tval[inter + 1]);
    for (int i = inter + 2; i < NEQ - 1; i++) {
        ri = data->x[i - 1];
        riP = data->x[i - 1] + h;
        riM = data->x[i - 2];
        Tpval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], N2, N2) / Cp(N2) / rho(N2, Tval[i]);// - функция справа от производных
    }
    ri = data->x[NEQ - 2];
    riP = data->x[NEQ - 2] + h;
    riM = data->x[NEQ - 3];
    Tpval[NEQ - 1] = Func_therm(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], N2, N2) / Cp(N2) / rho(N2, Tval[NEQ - 1]);
    /* Integration limits */
    t0 = ZERO;
    tout1 = tau;

    PrintHeader(rtol, avtol, T);

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);
    retval = IDAInit(mem, resEVEPORATION, t0, T, Tp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    IDASetMaxStep(mem, tau);
    /* Call IDASVtolerances to set tolerances */
    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    data->mem = mem;

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(T, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine */
   // retval = IDASetJacFn(mem, jacrob);
    //if (check_retval(&retval, "IDASetJacFn", 1)) return(1);

    /* Create Newton SUNNonlinearSolver object. IDA uses a
     * Newton SUNNonlinearSolver by default, so it is unecessary
     * to create it and attach it. It is done in this example code
     * solely for demonstration purposes. */
    NLS = SUNNonlinSol_Newton(T, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */
    ri = data->x[inter];
    riP = data->x[inter + 1];
    riM = data->x[inter - 1];
    iout = 0; tout = tout1;
    data->tout = 0;
    std::ofstream outf;
    outf.open("dot_M.txt");
    std::ofstream outp;
    outp.open("r(t).txt");
    while (1) {
        //std::cout << "inter  " << data->x[inter - 1] + p * h << std::endl;
        p = data->p;
        inter = data->inter;
        ri = data->x[inter - 1] + p * h;
        retval = IDASolve(mem, tout, &tret, T, Tp, IDA_NORMAL);
        // std::cout << " t = " << tout << " dot_M = " << Tval[N_d] <<"p = " << p << std::endl;
        if (check_retval(&retval, "IDASolve", 1)) return(1);

        if (retval == IDA_ROOT_RETURN) {
            retvalr = IDAGetRootInfo(mem, rootsfound);
            check_retval(&retvalr, "IDAGetRootInfo", 1);
            PrintRootInfo(rootsfound[0], rootsfound[1]);
        }
        //outf << tout << "," << lambda(H2O) * (Tval[inter + 1] - Tval[inter])/(riP - riM) << std::endl;
        if (retval == IDA_SUCCESS) {
            ri = data->x[inter - 1] + p * h;
            outf << data->t << "," << Tval[N_d] << std::endl;
            outp << data->t << "," << ri * ri << std::endl;
            iout++;
            data->t = tout;
            tout += tau;
            //inter = data->inter;

        }
        double delta = TOUT - tout;
        if (delta < 1.e-7) break;
        //std::cout << "\n\n\n\n\n\n\n\n";
    }
    for (int i = 0; i < NEQ; i++) {
        inter_T[i] = Tval[i];
    }
    outf.close();
    outp.close();
    out_Inter(tout, Tval, data);
    out_T(tout, Tval, data);
    printf("\n\n");
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(T);
    N_VDestroy(Tp);
    SUNContext_Free(&ctx);

    return(retval);

}
int main(void)
{
    void* mem;
    realtype T[NEQ];
    realtype t0, tout1, tout, tret, h, ri, riM, riP, v;
    int iout, retval, retvalr;
    int inter = N_d - 1;
    int rootsfound[2];
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    realtype p = p0;
    SUNContext ctx;
    UserData data;
    std::vector<double> R_vec;
    v = V0; //скорость движения интерфейса
    h = L / (N - 1);
    mem = NULL;
    data = NULL;
    R_vec = mesh();
    // std::cout << "WIDTH + N_d - 2 = " << R_vec.at(WIDTH + N_d - 2) << "\n";
    data = (UserData)malloc(sizeof * data);
    if (data) {
        for (int i = 0; i < R_vec.size(); i++) {
            data->x[i] = R_vec.at(i);
        }
    }
    data->p = p;
    data->p_new = p0;
    data->step = false;
    data->inter = inter;
    data->t = 0.0;
    //inital_T(initT, data);
    heat_inter(T, data);
    evaporation(T, data);

    return 1;

}
/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */

 /*
  * Define the system residual function.
  */
  //Idagetcurenttime

int resHEAT(realtype tres, N_Vector T, N_Vector Tp, N_Vector rr, void* user_data)
{
    int inter;
    realtype* Tval, * Tpval, * rval;
    realtype h = L / (N - 1);
    realtype ri, riP, riM, p;
    realtype tprev;
    UserData data;
    SUNContext ctx;
    int retval;
    void* mem;
    data = (UserData)user_data;
    Tval = N_VGetArrayPointer(T);
    Tpval = N_VGetArrayPointer(Tp);
    rval = N_VGetArrayPointer(rr);
    mem = data->mem;
    /*std::cout << tres << "  tres" << std::endl;*/
    p = data->p;
    ri = data->x[0];
    riM = 0;
    riP = data->x[1];

    inter = data->inter;

    rval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[0]) * Tpval[0];
    for (int i = 1; i < inter; i++) {
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[i]) * Tpval[i];
    }

    /* p = data->p;*/
     /*std::cout << p << " = p" << std::endl;*/
    h = data->x[inter] - data->x[inter - 1];
    riP = data->x[inter - 1] + p * h;
    riM = data->x[inter - 1];
    ri = data->x[inter];
    rval[inter] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 2], Tval[inter - 1], Tval[inter], Tval[inter + 1], H2O, H2O, h, p) - Cp(H2O) * rho(H2O, Tval[inter]) * Tpval[inter];
    //std::cout << rval[inter] << "  inter - 1" << std::endl;
    /*rval[inter - 1] = Func_therm_inter(Tval[inter - 3], Tval[inter - 2], Tval[inter], Tval[inter + 2], Tval[inter + 3], H2O, N2, h, p);*/
    rval[inter + 1] = Func_therm_inter(Tval[inter - 2], Tval[inter - 1], Tval[inter + 1], Tval[inter + 3], Tval[inter + 4], H2O, N2, h, p);
    // std::cout << rval[inter + 1] <<"  " << p << std::endl;
    ri = data->x[inter + 1];
    riP = data->x[inter + 2];
    riM = data->x[inter - 1] + p * h;
    rval[inter + 2] = Func_therm_inter_P(riM, ri, riP, Tval[inter + 1], Tval[inter + 2], Tval[inter + 3], Tval[inter + 4], N2, N2, h, p) - Tpval[inter + 2] * Cp(N2) * rho(N2, Tval[inter + 2]);
    for (int i = inter + 3; i < NEQ - 1; i++) {
        //std::cout <<" i =  " << i << std::endl;
        ri = data->x[i - 1];
        riP = data->x[i - 1] + h;
        riM = data->x[i - 2];
        rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], N2, N2) - Cp(N2) * rho(N2, Tval[i]) * Tpval[i];
    }
    ri = data->x[N - 2];
    riP = data->x[N - 2] + h;
    riM = data->x[N - 3];

    rval[NEQ - 1] = Func_therm(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], N2, N2) - Cp(N2) * rho(N2, Tval[NEQ - 1]) * Tpval[NEQ - 1];
    return 0;
}

int resHEAT_inital(realtype tres, N_Vector T, N_Vector Tp, N_Vector rr, void* user_data)
{
    int inter;
    realtype* Tval, * Tpval, * rval;
    realtype h = L / (N - 1);
    realtype ri, riP, riM;
    UserData data;
    SUNContext ctx;
    int retval;
    void* mem;
    data = (UserData)user_data;
    Tval = N_VGetArrayPointer(T);
    Tpval = N_VGetArrayPointer(Tp);
    rval = N_VGetArrayPointer(rr);
    mem = data->mem;
    //std::cout << tres << "  tres" << std::endl;
    ri = data->x[0];
    riM = 0;
    riP = data->x[1];

    inter = data->inter;
    rval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[0]) * Tpval[0];
    for (int i = 1; i < inter; i++) {
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[i]) * Tpval[i];
    }
    h = data->x[inter] - data->x[inter - 1];
    riP = data->x[inter + 1];
    riM = data->x[inter - 1];
    ri = data->x[inter];
    rval[inter] = Func_therm(riM, ri, riP, Tval[inter - 1], Tval[inter], Tval[inter + 1], H2O, N2) - Cp(H2O) * rho(H2O, Tval[inter]) * Tpval[inter];
    ri = data->x[inter + 1];
    riP = data->x[inter + 2];
    riM = data->x[inter];
    rval[inter + 1] = Func_therm(riM, ri, riP, Tval[inter], Tval[inter + 1], Tval[inter + 2], N2, N2) - Cp(N2) * rho(N2, Tval[inter + 1]) * Tpval[inter + 1];
    for (int i = inter + 2; i < NEQ - 2; i++) {
        //std::cout <<" i =  " << i << std::endl;
        ri = data->x[i];
        riP = data->x[i + 1];
        riM = data->x[i - 1];
        rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], N2, N2) - Cp(N2) * rho(N2, Tval[i]) * Tpval[i];
    }
    ri = data->x[N - 2];
    riP = data->x[N - 2] + h;
    riM = data->x[N - 3];

    rval[NEQ - 2] = Func_therm(riM, ri, riP, Tval[NEQ - 3], Tval[NEQ - 2], Tval[NEQ - 2], N2, N2) - Cp(N2) * rho(N2, Tval[NEQ - 1]) * Tpval[NEQ - 2];
    return 0;
}

int resEVEPORATION(realtype tres, N_Vector T, N_Vector Tp, N_Vector rr, void* user_data) {
    int inter;
    realtype* Tval, * Tpval, * rval;
    realtype h = L / (N - 1);
    double dot_M;
    realtype ri, riP, riM, p;

    UserData data;
    SUNContext ctx;
    int retval;
    void* mem;
    data = (UserData)user_data;
    Tval = N_VGetArrayPointer(T);
    Tpval = N_VGetArrayPointer(Tp);
    rval = N_VGetArrayPointer(rr);
    mem = data->mem;
    inter = data->inter;
    p = data->p;
    h = data->x[inter] - data->x[inter - 1];
    ri = data->x[inter - 1] + p * h;
    dot_M = Tval[N_d];
    p = get_p(tres, dot_M, ri, h, data, mem);
    ri = data->x[0];
    riM = 0;
    riP = data->x[1];
    inter = data->inter;
    rval[0] = Func_therm(riM, ri, riP, Tval[0], Tval[0], Tval[1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[0]) * Tpval[0];
    if (inter == N_d - 1) {
        for (int i = 1; i < inter; i++) {
            ri = data->x[i];
            riP = data->x[i + 1];
            riM = data->x[i - 1];
            rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[i]) * Tpval[i];
            //std::cout<<"rval[i] = " << rval[i] <<"  " << i << std::endl;
        }
        h = data->x[inter] - data->x[inter - 1];
        ri = data->x[inter - 1] + p * h;
        dot_M = Tval[N_d];
        riP = data->x[inter - 1] + p * h;
        riM = data->x[inter - 1];
        ri = data->x[inter];
        rval[inter] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 2], Tval[inter - 1], Tval[inter], T_vep, H2O, H2O, h, p) - Cp(H2O) * rho(H2O, T_vep) * Tpval[inter];
        rval[N_d] = Func_therm_evoparation_inter(Tval[inter - 2], Tval[inter - 1], T_vep, Tval[inter + 3], Tval[inter + 4], dot_M, ri, H2O, N2, h, p);
        std::cout << " rval[Nd] = " << rval[N_d] << " p = " << p << std::endl;
        //std::cout << "p = " << p << std::endl;
        //std::cout << "Nd_rval = " << rval[N_d] << std::endl;
        //std::cout << " rval[Nd] = " << rval[N_d] << " p = " << p << std::endl;
        ri = data->x[inter + 1];
        riP = data->x[inter + 2];
        riM = data->x[inter - 1] + p * h;
        /* std::cout << "ri = " << riM << std::endl;*/
        rval[inter + 2] = Func_therm_inter_evoparation_P(riM, ri, riP, T_vep, Tval[inter + 2], Tval[inter + 3], Tval[inter + 4], dot_M, N2, N2, h, p) - Tpval[inter + 2] * Cp(N2) * rho(N2, Tval[inter + 2]);
        for (int i = inter + 3; i < NEQ - 1; i++) {
            if (i != N_d) {
                ri = data->x[i - 1];
                riP = data->x[i - 1] + h;
                riM = data->x[i - 2];
                rval[i] = Func_therm_evoparation(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[i]) * Tpval[i];
            }
        }
        ri = data->x[N - 2];
        riP = data->x[N - 2] + h;
        riM = data->x[N - 3];
        rval[NEQ - 1] = Func_therm_evoparation(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[NEQ - 1]) * Tpval[NEQ - 1];

        return 0;
    }
    if (inter == N_d - 2) {
        for (int i = 1; i < inter; i++) {
            ri = data->x[i];
            riP = data->x[i + 1];
            riM = data->x[i - 1];
            rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[i]) * Tpval[i];
            //std::cout<<"rval[i] = " << rval[i] <<"  " << i << std::endl;
        }

        /* p = data->p;*/
        // std::cout << p << " = p" << std::endl;
        h = data->x[inter] - data->x[inter - 1];
        ri = data->x[inter - 1] + p * h;
        dot_M = Tval[N_d];
        riP = data->x[inter - 1] + p * h;
        riM = data->x[inter - 1];
        ri = data->x[inter];
        rval[inter - 1] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 3], Tval[inter - 2], Tval[inter - 1], T_vep, H2O, H2O, h, p) - Cp(H2O) * rho(H2O, Tval[inter - 1]) * Tpval[inter - 1];
        riM = data->x[inter - 1] + p * h;
        ri = data->x[inter + 1];
        riP = data->x[inter + 2];
        rval[inter] = Func_therm_inter_evoparation_P(riM, ri, riP, T_vep, Tval[inter], Tval[inter + 3], Tval[inter + 4], dot_M, N2, N2, h, p) - Tpval[inter] * Cp(N2) * rho(N2, Tval[inter]);
        ri = data->x[inter - 1] + p * h;
        rval[N_d] = Func_therm_evoparation_inter(Tval[inter - 3], Tval[inter - 2], T_vep, Tval[inter + 3], Tval[inter + 4], dot_M, ri, H2O, N2, h, p);
        std::cout << " rval[Nd] = " << rval[N_d] << " p = " << p << std::endl;
        rval[N_d + 1] = Func_therm_evoparation(riM, ri, riP, Tval[N_d - 1], Tval[N_d + 1], Tval[N_d + 2], dot_M, N2, N2);
        // std::cout << "dot_M = " << dot_M << std::endl;
        ri = data->x[inter + 1];
        riP = data->x[inter + 2];
        riM = data->x[inter - 1] + p * h;
        /* std::cout << "ri = " << riM << std::endl;*/
        for (int i = inter + 1; i < NEQ - 1; i++) {
            if ((i != N_d) && (i != (N_d + 1))) {
                ri = data->x[i - 1];
                riP = data->x[i - 1] + h;
                riM = data->x[i - 2];
                rval[i] = Func_therm_evoparation(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[i]) * Tpval[i];
            }
        }
        ri = data->x[N - 2];
        riP = data->x[N - 2] + h;
        riM = data->x[N - 3];
        rval[NEQ - 1] = Func_therm_evoparation(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[NEQ - 1]) * Tpval[NEQ - 1];

        return 0;
        if (inter < N_d - 2) {
            for (int i = 1; i < inter; i++) {
                ri = data->x[i];
                riP = data->x[i + 1];
                riM = data->x[i - 1];
                rval[i] = Func_therm(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], H2O, H2O) - Cp(H2O) * rho(H2O, Tval[i]) * Tpval[i];
            }
            h = data->x[inter] - data->x[inter - 1];
            ri = data->x[inter - 1] + p * h;
            dot_M = Tval[N_d];
            p = get_p(tres, dot_M, ri, h, data, mem);
            riP = data->x[inter - 1] + p * h;
            riM = data->x[inter - 1];
            ri = data->x[inter];
            rval[inter - 1] = Func_therm_inter_M(riM, ri, riP, Tval[inter - 3], Tval[inter - 2], Tval[inter - 1], T_vep, H2O, H2O, h, p) - Cp(H2O) * rho(H2O, Tval[inter - 1]) * Tpval[inter - 1];
            riM = data->x[inter - 1] + p * h;
            ri = data->x[inter + 1];
            riP = data->x[inter + 2];
            rval[inter] = Func_therm_inter_evoparation_P(riM, ri, riP, T_vep, Tval[inter], Tval[inter + 2], Tval[inter + 3], dot_M, N2, N2, h, p) - Tpval[inter] * Cp(N2) * rho(N2, Tval[inter]);
            ri = data->x[inter - 1] + p * h;
            rval[N_d - 1] = Func_therm_evoparation(riM, ri, riP, Tval[N_d - 2], Tval[N_d - 1], Tval[N_d + 1], dot_M, N2, N2);
            rval[N_d] = Func_therm_evoparation_inter(Tval[inter - 3], Tval[inter - 2], T_vep, Tval[inter + 2], Tval[inter + 3], dot_M, ri, H2O, N2, h, p);
            rval[N_d + 1] = Func_therm_evoparation(riM, ri, riP, Tval[N_d - 1], Tval[N_d + 1], Tval[N_d + 2], dot_M, N2, N2);
            ri = data->x[inter + 1];
            riP = data->x[inter + 2];
            riM = data->x[inter - 1] + p * h;
            for (int i = inter + 1; i < NEQ - 1; i++) {
                if ((i != N_d) && (i != (N_d + 1)) && (i != (N_d - 1))) {
                    ri = data->x[i - 1];
                    riP = data->x[i - 1] + h;
                    riM = data->x[i - 2];
                    rval[i] = Func_therm_evoparation(riM, ri, riP, Tval[i - 1], Tval[i], Tval[i + 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[i]) * Tpval[i];
                }
            }
            ri = data->x[N - 2];
            riP = data->x[N - 2] + h;
            riM = data->x[N - 3];
            rval[NEQ - 1] = Func_therm_evoparation(riM, ri, riP, Tval[NEQ - 2], Tval[NEQ - 1], Tval[NEQ - 1], dot_M, N2, N2) - Cp(N2) * rho(N2, Tval[NEQ - 1]) * Tpval[NEQ - 1];

            return 0;
        }
    }
}
/*
 * Root function routine. Compute functions g_i(t,y) for i = 0,1.
 */


 /*
  * Define the Jacobian function.
  */

  /*
   *--------------------------------------------------------------------
   * Private functions
   *--------------------------------------------------------------------
   */

   /*
    * Print first lines of output (problem description)
    */

static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)
{
    realtype* atval, * yval;

    atval = N_VGetArrayPointer(avtol);
    yval = N_VGetArrayPointer(y);

    printf("\nidaRoberts_dns: Robertson kinetics DAE serial example problem for IDA\n");
    printf("         Three equation chemical kinetics problem.\n\n");
    printf("Linear solver: DENSE, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
        yval[0], yval[1], yval[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%g %g %g)\n",
        yval[0], yval[1], yval[2]);
#else
    printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
        rtol, atval[0], atval[1], atval[2]);
    printf("Initial conditions y0 = (%g %g %g)\n",
        yval[0], yval[1], yval[2]);
#endif
    printf("Constraints and id not used.\n\n");
    printf("-----------------------------------------------------------------------\n");
    printf("  t             y1           y2           y3");
    printf("      | nst  k      h\n");
    printf("-----------------------------------------------------------------------\n");
}

/*
 * Print Output
 */
static void PrintOutput(void* mem, realtype t, N_Vector y)
{
    realtype* yval;
    int retval, kused;
    long int nst;
    realtype hused;

    yval = N_VGetArrayPointer(y);

    retval = IDAGetLastOrder(mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1);
    retval = IDAGetNumSteps(mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1);
    retval = IDAGetLastStep(mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
    printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n",
        t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

static void PrintRootInfo(int root_f1, int root_f2)
{
    printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);
    return;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void* returnvalue, const char* funcname, int opt)
{
    int* retval;
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }
    else if (opt == 1) {
        /* Check if retval < 0 */
        retval = (int*)returnvalue;
        if (*retval < 0) {
            fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *retval);
            return(1);
        }
    }
    else if (opt == 2 && returnvalue == NULL) {
        /* Check if function returned NULL pointer - no memory allocated */
        fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    return(0);
}

/* compare the solution at the final time 4e10s to a reference solution computed
   using a relative tolerance of 1e-8 and absoltue tolerance of 1e-14 */
static int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol)
{
    int      passfail = 0;        /* answer pass (0) or fail (1) retval */
    N_Vector ref;               /* reference solution vector        */
    N_Vector ewt;               /* error weight vector              */
    realtype err;               /* wrms error                       */

    /* create reference solution and error weight vectors */
    ref = N_VClone(y);
    ewt = N_VClone(y);

    /* set the reference solution data */
    NV_Ith_S(ref, 0) = RCONST(5.2083474251394888e-08);
    NV_Ith_S(ref, 1) = RCONST(2.0833390772616859e-13);
    NV_Ith_S(ref, 2) = RCONST(9.9999994791631752e-01);

    /* compute the error weight vector, loosen atol */
    N_VAbs(ref, ewt);
    N_VLinearSum(rtol, ewt, RCONST(10.0), atol, ewt);
    if (N_VMin(ewt) <= ZERO) {
        fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
        return(-1);
    }
    N_VInv(ewt, ewt);

    /* compute the solution error */
    N_VLinearSum(ONE, y, -ONE, ref, ref);
    err = N_VWrmsNorm(ref, ewt);

    /* is the solution within the tolerances? */
    passfail = (err < ONE) ? 0 : 1;

    if (passfail) {
        //fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%"GSYM"\n\n", err);
    }

    /* Free vectors */
    N_VDestroy(ref);
    N_VDestroy(ewt);

    return(passfail);
}

