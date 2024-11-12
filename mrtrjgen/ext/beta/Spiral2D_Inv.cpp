#include "Spiral2D_Inv.h"
#include <stdexcept>

Spiral2D_Inv::Spiral2D_Inv() {}

Spiral2D_Inv::~Spiral2D_Inv() {}

Spiral2D_Inv& Spiral2D_Inv::operator=(Spiral2D_Inv &o)
{
	return *this;
}

void Spiral2D_Inv::Update(std::vector<double> *pvdPara)
{
    double dNp, dU, dTht0, dKmax, dS, dDt, dOv;
    std::list<double*> lpdPtrPara(0);
    lpdPtrPara.push_back(&dNp);
    lpdPtrPara.push_back(&dU);
    lpdPtrPara.push_back(&dTht0);
    lpdPtrPara.push_back(&dKmax);
    lpdPtrPara.push_back(&dS);
    lpdPtrPara.push_back(&dDt);
    lpdPtrPara.push_back(&dOv);
    std::vector<double*> vpdPtrPara = std::vector<double*>(lpdPtrPara.begin(), lpdPtrPara.end());
    ExtractPara(pvdPara, &vpdPtrPara);

    std::list<double> ldKx(0), ldKy(0);
    double dDtOv = dDt/dOv;
#define PUSHBACK_KXY() \
    {\
        ldKx.push_back(dRho*cos(dTht + dTht0));\
        ldKy.push_back(dRho*sin(dTht + dTht0));\
    }

    // k_previous
    double dKx_Prev1, dKx_Prev2, dKy_Prev1, dKy_Prev2;

    // k[0]
    double dRho = dKmax;
    double dTht = 0;
    dKx_Prev2 = dRho*cos(dTht + dTht0);
    dKy_Prev2 = dRho*sin(dTht + dTht0);
    PUSHBACK_KXY();

    // k[1]
    dTht += dS*dDtOv*dDtOv/dKmax;
    dRho += -(0.5e0)/(2e0*M_PI)*(dU)/(dNp/2e0)*dTht;
    dKx_Prev1 = dRho*cos(dTht + dTht0);
    dKy_Prev1 = dRho*sin(dTht + dTht0);
    PUSHBACK_KXY();

    // k[2] ~ k[inf]
    double dD1Tht = dTht/dDtOv, dD2Tht = 0;
    while (dRho > 0)
    {
        double dKx0, dKy0, dS0;
        {
            dD2Tht = SovD2Tht(dTht, dD1Tht, dNp, dU, dS, 1e0);
            double dD1Tht_Temp = dD1Tht + dD2Tht*dDtOv;
            double dTht_Temp = dTht + dD1Tht_Temp*dDtOv;
            double dRho_Temp = dRho - (0.5e0)/(2e0*M_PI)*(dU)/(dNp/2e0)*dD1Tht*dDtOv;
            dKx0 = dRho_Temp*cos(dTht_Temp);
            dKy0 = dRho_Temp*sin(dTht_Temp);
            double dSx = (dKx0 - 2*dKx_Prev1 + dKx_Prev2)/(dDtOv*dDtOv);
            double dSy = (dKy0 - 2*dKy_Prev1 + dKy_Prev2)/(dDtOv*dDtOv);
            dS0 = sqrt(dSx*dSx + dSy*dSy);
        }
        double dKx1, dKy1, dS1;
        {
            dD2Tht = SovD2Tht(dTht, dD1Tht, dNp, dU, dS, -1e0);
            double dD1Tht_Temp = dD1Tht + dD2Tht*dDtOv;
            double dTht_Temp = dTht + dD1Tht_Temp*dDtOv;
            double dRho_Temp = dRho - (0.5e0)/(2e0*M_PI)*(dU)/(dNp/2e0)*dD1Tht*dDtOv;
            dKx1 = dRho_Temp*cos(dTht_Temp);
            dKy1 = dRho_Temp*sin(dTht_Temp);
            double dSx = (dKx1 - 2*dKx_Prev1 + dKx_Prev2)/(dDtOv*dDtOv);
            double dSy = (dKy1 - 2*dKy_Prev1 + dKy_Prev2)/(dDtOv*dDtOv);
            dS1 = sqrt(dSx*dSx + dSy*dSy);
        }

        if (fabs(dS0 - dS) <= fabs(dS1 - dS) || ldKx.size() < 20)
        {
            printf("flag0\n");
            dD2Tht = SovD2Tht(dTht, dD1Tht, dNp, dU, dS, 1e0);
        }
        else
        {
            printf("flag1\n");
            dD2Tht = SovD2Tht(dTht, dD1Tht, dNp, dU, dS, -1e0);
        }
        dD1Tht += dD2Tht*dDtOv;
        dTht += dD1Tht*dDtOv;
        dRho -= (0.5e0)/(2e0*M_PI)*(dU)/(dNp/2e0)*dD1Tht*dDtOv;

        printf("dTht: %lf\n", dTht);

        dKx_Prev2 = dKx_Prev1;
        dKy_Prev2 = dKy_Prev1;
        dKx_Prev1 = dRho*cos(dTht + dTht0);
        dKy_Prev1 = dRho*sin(dTht + dTht0);
        PUSHBACK_KXY();
    }
#undef PUSHBACK_KXY
    m_lNpt = ldKx.size();
    
    // copy to member variables
    m_vdKx.assign(ldKx.begin(), ldKx.end());
    m_vdKy.assign(ldKy.begin(), ldKy.end());

    // update gradient vector
    m_vdGx.resize(m_lNpt); m_vdGy.resize(m_lNpt);
    m_vdGx[0] = 0e0; m_vdGy[0] = 0e0;
    for (int64_t lIdxPt = 1; lIdxPt < m_lNpt; ++lIdxPt)
    {
        m_vdGx[lIdxPt] = (m_vdKx[lIdxPt] - m_vdKx[lIdxPt-1])/dDt;
        m_vdGy[lIdxPt] = (m_vdKy[lIdxPt] - m_vdKy[lIdxPt-1])/dDt;
    }
}

double Spiral2D_Inv::SovD2Tht(double dD0Tht, double dD1Tht, double dNp, double dU, double dS, double dSign)
{
    double dA = (-0.5*M_PI*dD0Tht*dNp*dU + 0.25*pow(M_PI, 2)*pow(dNp, 2) + 0.25*pow(dU, 2)*(pow(dD0Tht, 2) + 1))/(pow(M_PI, 2)*pow(dNp, 2));
    double dB = 0.5*pow(dD1Tht, 2)*dU*(dD0Tht*dU - M_PI*dNp)/(pow(M_PI, 2)*pow(dNp, 2));
    double dC = (-0.5*M_PI*dD0Tht*pow(dD1Tht, 4)*dNp*dU + pow(dD1Tht, 4)*pow(dU, 2)*(0.25*pow(dD0Tht, 2) + 1.0) + pow(M_PI, 2)*pow(dNp, 2)*(0.25*pow(dD1Tht, 4) - 1.0*pow(dS, 2)))/(pow(M_PI, 2)*pow(dNp, 2));

    return SovQuadEq(dA, dB, dC, dSign);
}

double Spiral2D_Inv::SovD1Tht(double dD0Tht, double dNp, double dU, double dS)
{
    double dA = (-0.5*M_PI*dD0Tht*dNp*dU + 0.25*pow(M_PI, 2)*pow(dNp, 2) + pow(dU, 2)*(0.25*pow(dD0Tht, 2) + 1.0))/(pow(M_PI, 2)*pow(dNp, 2));
    double dB = 0;
    double dC = -pow(dS, 2);

    return sqrt(fabs(SovQuadEq(dA, dB, dC)));
}
