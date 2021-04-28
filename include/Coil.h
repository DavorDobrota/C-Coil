#ifndef GENERAL_COIL_PROGRAM_COIL_H
#define GENERAL_COIL_PROGRAM_COIL_H

#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <windows.h>
#include <chrono>
#include <functional>

#include "../include/ctpl.h"
#include "../include/hardware_acceleration.h"

using namespace std;

#define Type double

//constants
const double Pi = 3.14159265357989323;
const double Mi = 0.0000001;

const Type Moderator = 5695.3125;

class Coil;

class thread_pool : public ctpl::thread_pool
{
public:
    thread_pool(){Init();}

    thread_pool(int n){Init(n);}

    void Init(int n = 0)
    {
        if(!n)
        {
            SYSTEM_INFO info;
            GetSystemInfo(&info);
            resize(info.dwNumberOfProcessors);
        }
        else
            resize(n);
    }

    void Set_priority(int priority)
    {
        for(int it = 0; it < this->size(); it++)
            for(int it = 0; it < this->size(); it++)
                SetThreadPriority(reinterpret_cast<HANDLE>(this->get_thread(it).native_handle()), priority);
    }
};

struct thread_param
{
    Coil *c;
    Type r;
    Type Theta;
    Type dR_start;
    Type dR_stop;
    Type res;
    bool completed;
};


class Coil {

public:

    Type j;		//Current density, calculated from current and d
    Type I;		//Current in Amperes
    Type R;		//Internal coio radius
    Type a;		//Coil thickness
    Type b;		//Coil length
    Type d;		//Coil winding wire thickness (approximate)
    int nA;		//Number of precision increments over a
    Type IncA;	//Precision increment over a
    int nB;		//Number of precision increments over b
    Type IncB;	//Precision increment over b
    int nFi;	//Number of precision increments over angle Phi
    Type IncFi;	//Precision increment over angle, usually fixed to 12
    int N;		//Number of windings of a coil
    bool isEM;	//Determines if the coil should generate an electric field
    Type f;		//Frequency of sinusoidal oscillation in EM coils
    Type Ro;	//Resistivity of a coil, material dependent (1.63e-8 for copper)
    Type L;		//Self inductance of a coil, calculated later
    Type Res;	//Resistance, real part of impedence, skin effect compensated for AC fields
    Type Reac;	//Reactance, imaginary part of impedence, dependent on f and L
    Type Imp;	//Impedence, total dissipative potential in AC fields
    Type mM;	//Coil magnetic moment

    Type DWeights[4][4];	//static weights matrix used for 2D Boole integral smoothening
    vector<Type> FiWeights;	//contains sequence of weights used for 1D Boole smoothening
    vector<Type> EdgeList1; //contains edge weight used for 2D Boole smoothening
    vector<Type> EdgeList2;	//contains edge weight used for 2D Boole smoothening

    //Constructor with all coil characteristics accounted for
    Coil (Type _I, Type _R, Type _a, Type _b, Type _d, int _nA, int _nB, int _nFi, int _N, bool _isEM, Type _f, Type _Ro){

        j = _I*_N/(_a*_b);
        I = _I;
        R = _R;
        a = _a;
        b = _b;
        d = _d;
        nA = _nA;
        IncA = a / nA;
        nB = _nB;
        IncB = b / nB;
        nFi = _nFi;
        IncFi = Pi / _nFi;
        N = _N;
        isEM = _isEM;
        Ro = _Ro;

        if (_isEM==true){
            f = _f;
        }
        else{
            f = 0;
        }

        //initialising the static weights matrix
        for (int ii = 0; ii < 4; ++ii){
            for (int ij = 0; ij < 4; ++ij){
                if (ij%4 == 0 && ii%4 == 0){
                    DWeights[ii][ij] = 196 / Moderator;
                }
                else if (ii%4 == 0 && ij%4 == 1 || ii%4 == 1 && ij%4 == 0 || ii%4 == 3 && ij%4 == 0 || ii%4 == 0 && ij%4 == 3){
                    DWeights[ii][ij] = 448 / Moderator;
                }
                else if (ii%4 == 2 && ij%4 == 0 || ii%4 == 0 && ij%4 == 2){
                    DWeights[ii][ij] = 168 / Moderator;
                }
                else if (ii%4 == 1 && ij%4 == 1 || ii%4 == 3 && ij%4 == 3 || ii%4 == 3 && ij%4 == 1 || ii%4 == 1 && ij%4 == 3){
                    DWeights[ii][ij] = 1024 / Moderator;
                }
                else if (ii%4 == 2 && ij%4 == 1 || ii%4 == 1 && ij%4 == 2 || ii%4 == 3 && ij%4 == 2 || ii%4 == 2 && ij%4 == 3){
                    DWeights[ii][ij] = 384 / Moderator;
                }
                else if (ij%4 == 2 && ii%4 == 2){
                    DWeights[ii][ij] = 144 / Moderator;
                }
            }
        }

        initialiseLists();

        /*	SelfInductanceCalc();
            ResistanceCalc();
            ReactanceCalc();
            ImpedenceCalc();
            mMCalc();*/
    }

    //Default constructor
    Coil(){}

    void initialiseLists(){

        FiWeights.resize(0);
        EdgeList1.resize(0);
        EdgeList2.resize(0);

        for (int k=0; k<=nFi; k++){
            if (k==nFi || k==0){		FiWeights.push_back(7);}
            else if (k%4==0){			FiWeights.push_back(14);}
            else if (k%4==1 || k%4==3){	FiWeights.push_back(32);}
            else if (k%4==2){			FiWeights.push_back(12);}
        }

        for (int m=0; m<=nA; m++){
            if (m==0 || m==nA){		EdgeList1.push_back(49/Moderator);}
            else if (m%2 == 1){		EdgeList1.push_back(224/Moderator);}
            else if (m%4 == 0){		EdgeList1.push_back(98/Moderator);}
            else if (m%4 == 2){		EdgeList1.push_back(84/Moderator);}
        }

        for (int n=0; n<=nB; n++){
            if (n==0 || n==nB){		EdgeList2.push_back(49/Moderator);}
            else if (n%2 == 1){		EdgeList2.push_back(224/Moderator);}
            else if (n%4 == 0){		EdgeList2.push_back(98/Moderator);}
            else if (n%4 == 2){		EdgeList2.push_back(84/Moderator);}
        }
        return;
    }

    void newCurrent(Type _I){
        I = _I;
        j = _I * N / (a*b);
        mMCalc();
        return;
    }

    void setIncrements(int _nA, int _nB, int _nFi){
        nA = _nA;
        IncA = a / _nA;
        nB = _nB;
        IncB = b / _nB;
        nFi = _nFi;
        IncFi = Pi / _nFi;

        initialiseLists();
        return;
    }

    static long unsigned int BhAsyncThread(int id, void *in)
    {
        thread_param &p = *(thread_param*)in;
        const Coil &c = *p.c;

        Type FieldH = 0;
        Type G = 0, Dh = 0, E = 0, F = 0, eR = 0, er = 0;
        int ii = 0, ji = 0, ki = 0;
        Type Temp = 0;

        Type z = p.r * cos(p.Theta);
        Type q = p.r * sin(p.Theta);
        Type C = c.IncA * c.IncB * c.IncFi * c.j * Mi;
        Type bh = c.b/2;

        for (Type dR = p.dR_start; dR < p.dR_stop; dR += c.IncA){

            eR = c.R + dR;
            ii = int(round(dR/c.IncA));
            F = 2*eR * q;

            for (Type dr=-bh; dr<=bh; dr+=c.IncB){

                er = z + dr;
                ji = int(round((dr+bh)/c.IncB));
                Dh = C * eR * er;
                E = q*q + eR*eR + er*er;

                if (ii != c.nA && ii != 0 && ji != c.nB && ji != 0){
                    Temp = c.DWeights[ii%4][ji%4];
                }
                else if (ii == c.nA || ii == 0){
                    Temp = c.EdgeList2[ji];
                }
                else{
                    Temp = c.EdgeList1[ii];
                }

                for (Type dFi = 0; dFi <= Pi + c.IncFi/2; dFi += c.IncFi){

                    G = cos(dFi);
                    ki = int(round(dFi/c.IncFi));
                    FieldH += Temp * c.FiWeights[ki] * (Dh*G)/pow((E - F*G), 1.5);
                }
            }
        }
        p.res = FieldH;
        p.completed = true;
        return 0;
    }

    static long unsigned int BzAsyncThread(int id, void *in)
    {
        thread_param &p = *(thread_param*)in;
        const Coil &c = *p.c;

        Type FieldZ = 0;
        Type G = 0, Dz = 0, Cz = 0, E = 0, F = 0, eR = 0, er = 0;
        int ii = 0, ji = 0, ki = 0;
        Type Temp = 0;

        Type z = p.r * cos(p.Theta);
        Type q = p.r * sin(p.Theta);
        Type C = c.IncA * c.IncB * c.IncFi * c.j * Mi;
        Type bh=c.b/2;

        for (Type dR = p.dR_start; dR < p.dR_stop; dR += c.IncA){

            eR = c.R + dR;
            ii = int(round(dR/c.IncA));
            Cz = eR*eR;
            Dz = eR * q;
            F = 2*Dz;

            for (Type dr = -bh; dr <= bh; dr += c.IncB){

                er = z + dr;
                ji = int(round((dr+bh)/c.IncB));
                E = q*q + eR*eR + er*er;

                if (ii != c.nA && ii != 0 && ji != c.nB && ji != 0){
                    Temp = c.DWeights[ii%4][ji%4];
                }
                else if (ii == c.nA || ii == 0){
                    Temp = c.EdgeList2[ji];
                }
                else{
                    Temp = c.EdgeList1[ii];
                }

                for (Type dFi = 0; dFi <= Pi; dFi += c.IncFi){

                    G = cos(dFi);
                    ki=int(round(dFi/c.IncFi));
                    FieldZ += Temp * c.FiWeights[ki] * C * ((Cz - Dz*G)/(pow((E - F*G), 1.5)));
                }
            }
        }
        p.res = FieldZ;
        p.completed = true;
        return 0;
    }

    static long unsigned int EAsyncThread(int id, void *in)
    {
        thread_param &p = *(thread_param*)in;
        const Coil &c = *p.c;

        Type FieldE = 0;
        Type G = 0, E = 0, F = 0, eR = 0, er = 0, C = 0;
        int ii = 0, ji = 0, ki = 0;
        Type Temp = 0;

        Type z = p.r * cos(p.Theta);
        Type q = p.r * sin(p.Theta);
        Type C0 = c.IncA * c.IncB * c.IncFi * c.j * Mi;
        Type bh = c.b/2;

        for (Type dR = p.dR_start; dR < p.dR_stop; dR += c.IncA){

            eR = c.R + dR;
            C = C0 * eR;
            ii = int(round(dR/c.IncA));
            F = 2*eR * q;

            for (Type dr = -bh; dr <= bh; dr += c.IncB){

                er = z + dr;
                ji = int(round((dr+bh)/c.IncB));
                E = q*q + eR*eR + er*er;

                if (ii != c.nA && ii != 0 && ji != c.nB && ji != 0){
                    Temp = c.DWeights[ii%4][ji%4];
                }
                else if (ii == c.nA || ii == 0){
                    Temp = c.EdgeList2[ji];
                }
                else{
                    Temp = c.EdgeList1[ii];
                }

                for (Type dFi = 0; dFi <= Pi; dFi += c.IncFi){

                    G = cos(dFi);
                    ki = int(round(dFi/c.IncFi));
                    FieldE += Temp * c.FiWeights[ki] * C * ((G)/(sqrt(E - F*G)));
                }
            }
        }

        p.res = FieldE * 2*Pi * c.f;
        p.completed = true;
        return 0;
    }

    //Magnetic field calculation methods
    Type BhCalc(Type, Type);

    Type BxCalc(Type r, Type Theta, Type Alpha){

        return BhCalc(r, Theta) * cos(Alpha);
    }

    Type ByCalc(Type r, Type Theta, Type Alpha){

        return BhCalc(r, Theta) * sin(Alpha);
    }

    Type BzCalc(Type, Type);

    Type BCalc(Type r, Type Theta){

        Type Bh = BhCalc(r, Theta);
        Type Bz = BzCalc(r, Theta);

        return sqrt(Bh*Bh + Bz*Bz);
    }

    //Electric field calculation methods
    Type ECalc (Type, Type);

    Type ExCalc(Type r, Type Theta, Type Alpha){

        return -ECalc(r, Theta) * sin(Alpha);
    }

    Type EyCalc(Type r, Type Theta, Type Alpha){

        return ECalc(r, Theta) * cos(Alpha);
    }

    //Gradient field calculation methods
    Type GzCalc(Type r, Type Theta){

        Type dr = 0.0000001;
        Type rh = r * sin(Theta);
        Type rz = r * cos(Theta);

        Type r1 = sqrt(rh*rh + (rz - dr)*(rz - dr));
        Type r2 = sqrt(rh*rh + (rz + dr)*(rz + dr));
        Type Theta1 = acos((rz - dr) / r1);
        Type Theta2 = acos((rz + dr) / r2);

        Type B1 = BzCalc(r1, Theta1);
        Type B2 = BzCalc(r2, Theta2);

        return (B2 - B1)/(2*dr);
    }
    Type GhCalc(Type r, Type Theta){

        Type dr = 0.0000001;
        Type rh = r * sin(Theta);
        Type rz = r * cos(Theta);

        Type r1 = sqrt((rh - dr)*(rh - dr) + rz*rz);
        Type r2 = sqrt((rh + dr)*(rh + dr) + rz*rz);
        Type Theta1 = acos(rz / r1);
        Type Theta2 = acos(rz / r2);

        Type B1 = BhCalc(r1, Theta1);
        Type B2 = BhCalc(r2, Theta2);

        return (B2 - B1)/(2*dr);
    }

    Type GxCalc(Type r, Type Theta, Type Alpha){

        return GhCalc(r, Theta) * cos(Alpha);
    }
    Type GyCalc(Type r, Type Theta, Type Alpha){

        return GhCalc(r, Theta) * sin(Alpha);
    }

    //Inductance calculation, executed upon coil creation
    void SelfInductanceCalc(){

        for (Type id = 0; id <= b/2; id += d){

            for (Type jd = R; jd < R+a; jd += d){

                Type TempL = 0;

                Type r = sqrt((b/2 - id)*(b/2 - id) + jd*jd);
                Type Theta = asin(jd / r);

                Type G = 0, E = 0, F = 0, eR = 0, er = 0, C = 0;
                int ii = 0, ji = 0, ki = 0;
                Type Temp = 0;

                Type z = r * cos(Theta);
                Type q = r * sin(Theta);
                Type C0 = 2 * IncA * IncB * IncFi * Mi / (d*d);
                Type bh = b/2;

                for (Type dR = 0; dR <= a; dR += IncA){

                    eR = R + dR;
                    C = C0 * eR;
                    ii = int(round(dR/IncA));
                    F = 2*eR * q;

                    for (Type dr = -bh; dr <= bh; dr += IncB){

                        er = z + dr;
                        ji = int(round((dr+bh)/IncB));
                        E = q*q + eR*eR + er*er;

                        if (int((id)/IncB) == ji && int((jd-R)/IncA) == ii){

                            continue;
                        }
                        else{
                            if (ii != nA && ii != 0 && ji != nB && ji != 0){
                                Temp = DWeights[ii%4][ji%4];
                            }
                            else if (ii == nA || ii == 0){
                                Temp = EdgeList2[ji];
                            }
                            else{
                                Temp = EdgeList1[ii];
                            }
                            printf("%.10e\n", Temp*C);
                            for (Type dFi = 0; dFi <= Pi; dFi += IncFi){

                                ki = int(round(dFi/IncFi));
                                G = cos(dFi);
                                TempL += Temp * FiWeights[ki] * C * ((G)/(sqrt(E - F*G)));
                                printf("%.15f\n", (G)/(sqrt(E - F*G)));
                            }

                        }
                    }
                }
                //	printf("%2d %2d %2d %2d : %.15e %.15e %.15e\n", id, jd, ii, ji, C, E, F);
                if (isinf(TempL) == false){
                    L += TempL*2*Pi*jd;
                }
                else{
                    TempL = 0;
                }
            }
        }
        return;
    }

    void newSelfInductanceCalc(){

        L = 0;

        double TempL = 0;
        double rh, rz;
        double  E, F, G, cRh, cRz;

        double C = IncA * IncB * IncFi * Mi * N/(a*b) * (double) N/(nA*nB);
        double weight, temp;

        for (int id = 0; id <= nB; ++id){

            for (int jd = 0; jd <= nA; ++jd){

                TempL = 0;

                rh = jd * IncA + R;
                rz = id * IncB - b/2;

                if (jd != nA && jd != 0 && id != nB && id != 0){
                    weight = DWeights[jd%4][id%4];
                }
                else if (jd == nA || jd == 0){
                    weight = EdgeList2[id];
                }
                else{
                    weight = EdgeList1[jd];
                }

                for (int ia = 0; ia <= nA; ++ia){

                    for (int ib = 0; ib <= nB; ++ib){

                        cRh = IncA * ia + R;
                        cRz = IncB * ib + rz - b/2;

                        E = rh*rh + cRh*cRh + cRz*cRz;
                        F = 2 * cRh * rh;
                        temp = cRh;

                        if (id == ib && jd == ia){
                            continue;
                        }
                        else{

                            if (ia != nA && ia != 0 && ib != nB && ib != 0){
                                temp *= DWeights[ia%4][ib%4];
                            }
                            else if (ia == nA || ia == 0){
                                temp *= EdgeList2[ib];
                            }
                            else{
                                temp *= EdgeList1[ia];
                            }

                            for (int fi = 0; fi <= nFi; ++fi){

                                G = cos(fi * IncFi);
                                TempL += temp * FiWeights[fi] * C * ((G)/(sqrt(E - F*G)));
                                if (isinf(TempL)){
                                    printf("%d %d %d : %.10e %.10e %.10e %.10e %.10e\n", ib, ia, fi, temp, FiWeights[fi], C, G, E - F*G);
                                    TempL = 0;
                                }
                            }
                        }
                    }
                }
                //	printf("%2d %2d : %.15e %.15e %.15e\n", id, jd, C, E, F);
                //	printf("%2d %2d %.30f\n", id, jd, TempL);
                if (isinf(TempL) == false){
                    L += TempL * (45.0/4.0) * 2*Pi * rh;
                    printf("%2d %2d : %.15e %.15e %.15e\n", id, jd, C, E, F);
                }
                else{
                    TempL = 0;
                }
            }
        }
        return;
    }

    //Impedance and resistance calculation, executed upon coil creation
    void ResistanceCalc(){

        Type Resistance = N * (R + a/2);
        Type skin = sqrt(Ro / (Pi*f * Mi));

        Type normalCrossection = (d/2)*(d/2)*Pi;
        Type acCrossection = 2*Pi * (skin*skin * (exp(-(d/2)/skin) - 1) + skin * (d/2));

        Res = Resistance * (normalCrossection/acCrossection);

        return;
    }
    void ReactanceCalc(){
        Reac = L * 2*Pi*f;
        return;
    }
    void ImpedenceCalc(){
        Imp = sqrt(Res*Res + Reac*Reac);
        return;
    }

    //magnetic moment calculation,
    void mMCalc(){
        mM = Pi * I * N * (R*R + R*a + a*a/3);
        return;
    }

    //Interaction between coils: force
    Type ForceCalcBetweenCoilsAmpZ(Type dist, Coil coil2){

        Type Force = 0, B = 0, distance = 0, theta = 0, hComp = 0, zComp = 0;

        for (Type dR = coil2.IncA/2; dR < coil2.a; dR += coil2.IncA){

            hComp = coil2.R + dR;

            for (Type dr = 0; dr < coil2.b; dr += coil2.IncB){

                zComp = dist + dr - coil2.b/2 + coil2.IncB/2;
                distance = sqrt(hComp*hComp + zComp*zComp);
                theta = asin(hComp / distance);

                B = BhCalc(distance, theta);
                Force += coil2.IncA * coil2.IncB * j * 2*Pi * hComp * B;
            }
        }
        return Force;
    }
    Type ForceCalcBetweenCoilsGradZ(Type dist, Coil coil2){
        return -GzCalc(dist, 0.0) * coil2.mM;
    }
    Type ForceCalcBetweenCoilsMomentZ(Type dist, Coil coil2){
        return ((6*Mi * mM * coil2.mM) / pow(dist, 4));
    }

    //Interaction between coils: induction
    Type MutualInductanceCalc(Type dist, Coil sec, bool automatic, Type precision){

        //determining precision
        if (automatic == true){

            int minNF1, minNC1, minNC2;
            int m = 0, n = 0, o = 0;
            double nInc = 0.0, tempN = 0.0;

            if (precision >= 1.0 && precision <= 9.0){
                nInc = pow(2, 19 + precision);
                minNF1 = 12 + 4*int(floor(precision/2));
                minNC1 = 12 + 4*int(floor(precision/3));
                minNC2 = 12 + 4*int(floor(precision/4));
            }
            else if (precision > 9.0){
                nInc = pow(2, 28);
                minNF1 = 24;
                minNC1 = 24;
                minNC2 = 20;
            }
            else {
                nInc = pow(2, 20);
                minNF1 = 12;
                minNC1 = 12;
                minNC2 = 12;
            }

            int nF1, nC1, nC2;
            double f1Len, c1Len, c2Len;

            while (tempN < nInc){

                nF1 = minNF1 + 4*m;
                nC1 = minNC1 + 4*n;
                nC2 = minNC2 + 4*o;

                f1Len = Pi * (R+a/2) / nF1;
                c1Len = sqrt(a*b) / nC1;
                c2Len = sqrt(sec.a*sec.b) / nC2;

                if ((c1Len*c2Len) / (f1Len) <= 1){
                    m++;
                }
                else {

                    if (c2Len/c1Len <= 1){
                        n++;
                    } else{
                        o++;
                    }
                }
                tempN = nF1 * nC1*nC1 * nC2*nC2;
            }
            printf("%.2f %d %d %d %d %d\n", precision, int(nInc), nF1, nC1, nC2, int(tempN));

            setIncrements(nC1, nC1, nF1);
            sec.setIncrements(nC2, nC2, nF1);
        }

        double Mout = 0;
        int numFilaments = (sec.nA+1)*(sec.nB+1);

        vector <float> displacementArr;
        vector <float> thetaArr;
        vector <float> potentialArr(numFilaments);
        vector <double> dRArr;

        double dispH, dispZ, disp, temp;

        for (int ib = 0; ib <= sec.nB; ++ib){

            for (int ia = 0; ia <= sec.nA; ++ia){

                dispH = sec.R + ia*sec.IncA;
                dispZ = dist + ib*sec.IncB - sec.b/2;
                disp = sqrt(dispH*dispH + dispZ*dispZ);

                displacementArr.push_back(disp);
                thetaArr.push_back(acos(dispZ / disp));

                if (ib != sec.nB && ib != 0 && ia != sec.nA && ia != 0){
                    temp = sec.DWeights[ib%4][ia%4];
                }
                else if (ib == sec.nB || ib == 0){
                    temp = sec.EdgeList1[ia];
                }
                else{
                    temp = sec.EdgeList2[ib];
                }
                dRArr.push_back(dispH * temp);
            }
        }
        Calculate_hardware_accelerated_a(numFilaments, &thetaArr[0], &displacementArr[0], j, R, b, a, IncA, IncB, IncFi, nullptr, nullptr, &potentialArr[0]);

        for (int i = 0; i < numFilaments; ++i){
            Mout += potentialArr[i] * dRArr[i];
        }
        Mout = Mout * (45.0/4.0) / ((double)(sec.nA*sec.nB) / sec.N);

        return Mout;
    }

    Type MutualInductanceGeneralCalc(Coil sec, Type distance, Type offAxis, Type alpha, Type beta, bool automatic, Type precision){
        if (automatic == true){
            //setting precision
            int minNF1, minNF2, minNC1, minNC2;

            int m = 0, n = 0, o = 0, p = 0;
            double nInc = 0.0, tempN = 0.0;

            if (precision >= 1.0 && precision <= 9.0){
                nInc = pow(2, 23 + precision);
                minNF1 = 12 + 4*int(floor(precision/2));
                minNF2 = 12 + 4*int(floor(precision/2));
                minNC1 = 12 + 4*int(floor(precision/4));
                minNC2 = 12 + 4*int(floor(precision/4.5));
            }
            else if (precision > 9.0){
                nInc = pow(2, 32);
                minNF1 = 24;
                minNF2 = 24;
                minNC1 = 20;
                minNC2 = 20;
            }
            else {
                nInc = pow(2, 24);
                minNF1 = 12;
                minNF2 = 12;
                minNC1 = 12;
                minNC2 = 12;
            }

            int nF1, nF2, nC1, nC2;
            double f1Len, f2Len, c1Len, c2Len;

            while (tempN < nInc){

                nF1 = minNF1 + 4*m;
                nF2 = minNF2 + 4*n;
                nC1 = minNC1 + 4*o;
                nC2 = minNC2 + 4*p;

                f1Len = Pi * (R+a/2) / nF1;
                f2Len = Pi * (sec.R+sec.a/2) / nF2;
                c1Len = sqrt(a*b) / nC1;
                c2Len = sqrt(sec.a*sec.b) / nC2;

                if ((c1Len*c2Len) / (f1Len*f2Len) <= 1){
                    if (f2Len/f1Len <= 1){
                        m++;
                    } else{
                        n++;
                    }
                }
                else {
                    if (c2Len/c1Len <= 1){
                        o++;
                    } else{
                        p++;
                    }
                }
                tempN = nF1 * nF2 * nC1*nC1 * nC2*nC2;
            }

            //	printf("p = %.2f: %d %d %d %d %d %d\n", precision, int(nInc), nF1, nF2, nC1, nC2, int(tempN));
            setIncrements(nC1, nC1, nF1);
            sec.setIncrements(nC2, nC2, nF2);
        }

        double Mout = 0;
        int numFilaments = (sec.nA+1) * (sec.nB+1) * (sec.nFi+1);

        vector <float> displacementArr;
        vector <float> thetaArr;
        vector <float> potentialArr(numFilaments);
        vector <double> dRArr;

        vector <double> unitRingPointX, unitRingPointY, unitRingPointZ;

        double minFiRing, maxFiRing, incFiRing;

        if (beta == 0 || alpha == 0 || offAxis == 0){
            minFiRing = 0;
            maxFiRing = Pi;
            beta = 0;
        }
        else {
            minFiRing = 0;
            maxFiRing = 2*Pi;
        }
        incFiRing = (maxFiRing - minFiRing) / sec.nFi;

        for (double Fi = minFiRing; Fi <= maxFiRing; Fi += incFiRing){

            unitRingPointX.push_back(cos(beta)*cos(Fi) - sin(beta)*cos(alpha)*sin(Fi));
            unitRingPointY.push_back(sin(beta)*cos(Fi) + cos(beta)*cos(alpha)*sin(Fi));
            unitRingPointZ.push_back(sin(alpha)*sin(Fi));
        }

        double tempDispX, tempDispY, tempDispZ, tempDisp;
        double rho, weight;

        for (int ib = 0; ib <= sec.nB; ++ib){

            for (int ia = 0; ia <= sec.nA; ++ia){

                for (int ir = 0; ir <= sec.nFi; ++ir){

                    tempDispX = offAxis + (sec.R + ia*sec.IncA) * unitRingPointX[ir];
                    tempDispY = (sec.R + ia*sec.IncA) * unitRingPointY[ir];
                    tempDispZ = distance + ib*sec.IncB - sec.b/2 + (sec.R + ia*sec.IncA) * unitRingPointZ[ir];
                    tempDisp = sqrt(pow(tempDispX, 2) + pow(tempDispY, 2) + pow(tempDispZ, 2));

                    displacementArr.push_back(tempDisp);
                    thetaArr.push_back(acos(tempDispZ / tempDisp));

                    rho = atan2(tempDispY, tempDispX);

                    if (ib != sec.nB && ib != 0 && ia != sec.nA && ia != 0){
                        weight = sec.DWeights[ib%4][ia%4];
                    }
                    else if (ib == sec.nB || ib == 0){
                        weight = sec.EdgeList1[ia];
                    }
                    else{
                        weight = sec.EdgeList2[ib];
                    }
                    weight *= sec.FiWeights[ir];

                    dRArr.push_back(((sec.R + sec.IncA*ia) / (2*sec.nFi)) * weight * (cos(rho)*unitRingPointX[ir] + sin(rho)*unitRingPointY[ir]));
                }
            }
        }

        Calculate_hardware_accelerated_a(numFilaments, &thetaArr[0], &displacementArr[0], j, R, b, a, IncA, IncB, IncFi, nullptr, nullptr, &potentialArr[0]);

        for (int i = 0; i < numFilaments; ++i){
            Mout += (double)potentialArr[i] * dRArr[i];
        }
        Mout /= ((double)(sec.nA*sec.nB) / sec.N);

        return Mout;
    }

    Type VoltageOnSecondNoLoad (Type dist, Coil coil2){

        return MutualInductanceCalc(dist, coil2, true, 5.0)*2*Pi*f;
    }

    void setFreq(Type newFreq){
        f = newFreq;
        ResistanceCalc();
        ReactanceCalc();
        ImpedenceCalc();
        return;
    }

    //data coherence check
    bool operator!=(Coil &in)
    {
        if(in.j != j)
            return true;
        else if(in.I != I)
            return true;
        else if(in.R != R)
            return true;
        else if(in.a != a)
            return true;
        else if(in.b != b)
            return true;
        else if(in.d != d)
            return true;
        else if(in.nA != nA)
            return true;
        else if(in.IncA != IncA)
            return true;
        else if(in.nB != nB)
            return true;
        else if(in.IncB != IncB)
            return true;
        else if(in.nFi != nFi)
            return true;
        else if(in.IncFi != IncFi)
            return true;
        else if(in.N != N)
            return true;
        else if(in.isEM != isEM)
            return true;
        else if(in.f != f)
            return true;
        else if(in.Ro != Ro)
            return true;
        else if(in.FiWeights.size() != FiWeights.size())
            return true;
        else if(in.EdgeList1.size() != EdgeList1.size())
            return true;
        else if(in.EdgeList2.size() != EdgeList2.size())
            return true;

        return false;
    }

private:

    Type Calc(Type r, Type Theta, std::function<long unsigned int(int, void*)> func);
};

#endif //GENERAL_COIL_PROGRAM_COIL_H

