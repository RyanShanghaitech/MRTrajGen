#include "Spiral3D_A.h"

Spiral3D_A::Spiral3D_A() {}

Spiral3D_A::~Spiral3D_A() {}

Spiral3D_A& Spiral3D_A::operator=(Spiral3D_A &o)
{
	return *this;
}

double Spiral3D_A::SovD1Phi(double dD0Phi, double dNp, double dUTht, double dUPhi, double dTht0, double dPhi0, double dS)
{
    double dA = (1.0/32.0)*dUPhi*(16*pow(dD0Phi, 23.0/2.0)*dUTht*(dD0Phi*dUTht*pow(sin(dD0Phi + dPhi0), 2)*pow(sin(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + dD0Phi*dUTht*pow(sin(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)*pow(cos(dD0Phi + dPhi0), 2) + dUPhi*pow(sin(dD0Phi + dPhi0), 2)*pow(sin(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + 2*dUPhi*pow(sin(dD0Phi + dPhi0), 2)*pow(cos(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + dUPhi*pow(sin(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)*pow(cos(dD0Phi + dPhi0), 2) + 2*dUPhi*pow(cos(dD0Phi + dPhi0), 2)*pow(cos(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2)) + 4*pow(dD0Phi, 21.0/2.0)*pow(dUPhi, 2) + 24*pow(dD0Phi, 21.0/2.0)*pow(dUTht, 2)*pow(sin(M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + dTht0), 2) + 6*pow(dD0Phi, 19.0/2.0)*dUPhi*dUTht + pow(dD0Phi, 17.0/2.0)*pow(dUTht, 2) + 12*M_SQRT2*pow(dD0Phi, 11)*pow(dUTht, 2)*sqrt(dUPhi/dUTht)*sin(2*M_SQRT2*sqrt(dD0Phi)*sqrt(dUPhi/dUTht) + 2*dTht0))/(pow(M_PI, 2)*pow(dD0Phi, 23.0/2.0)*pow(dNp, 2)*dUTht);
    double dB = 0;
    double dC = -pow(dS, 2);

    return sqrt(fabs(SovQuadEq(dA, dB, dC)));
}

void Spiral3D_A::Update(std::vector<double> *pvdPara)
{
    double dNp, dUTht, dUPhi, dTht0, dPhi0, dKmax, dS, dDt;
    std::list<double*> ldParaDst(0);
    ldParaDst.push_back(&dNp);
    ldParaDst.push_back(&dUTht);
    ldParaDst.push_back(&dUPhi);
    ldParaDst.push_back(&dTht0);
    ldParaDst.push_back(&dPhi0);
    ldParaDst.push_back(&dKmax);
    ldParaDst.push_back(&dS);
    ldParaDst.push_back(&dDt);
    std::vector<double*> vdParaDst(ldParaDst.begin(), ldParaDst.end());
    ExtractPara(pvdPara, &vdParaDst);

    std::list<double> ldKx, ldKy, ldKz;
#define PUSHBACK_KXYZ() \
    {\
        ldKx.push_back(dRho*sin(dTht + dTht0)*cos(dPhi + dPhi0));\
        ldKy.push_back(dRho*sin(dTht + dTht0)*sin(dPhi + dPhi0));\
        ldKz.push_back(dRho*cos(dTht + dTht0));\
    }

    // k[0]
    double dRho = 0;
    double dTht = 0;
    double dPhi = 0;
    PUSHBACK_KXYZ();

    // k[1], since SovD1Phi() doesn't allow phi=0, we have to calculate k[1]
    dRho = (0.01*dS)*dDt*dDt/2;
    dTht = dRho*2*M_PI*dNp/dUTht;
    dPhi = (2*dRho*dRho*M_PI*M_PI*dNp*dNp)/(dUTht*dUPhi);
    PUSHBACK_KXYZ();

    // k[2] ~ k[inf]
    double dD1Phi = 0;
    while (dRho <= dKmax)
    {
        dD1Phi = SovD1Phi(dPhi, dNp, dUTht, dUPhi, dTht0, dPhi0, dS);
        dPhi += dD1Phi*dDt;

        dTht = sqrt((2*dUPhi)/(dUTht))*sqrt(dPhi);
        dRho = sqrt(dUTht*dUPhi/2)/(M_PI*dNp)*sqrt(dPhi);

        PUSHBACK_KXYZ();
    }
#undef PUSHBACK_KXYZ
    m_lNpt = ldKx.size();
    
    // copy to member variables
    m_vdKx.assign(ldKx.begin(), ldKx.end());
    m_vdKy.assign(ldKy.begin(), ldKy.end());
    m_vdKz.assign(ldKz.begin(), ldKz.end());

    // update gradient vector
    m_vdGx.resize(m_lNpt); m_vdGy.resize(m_lNpt); m_vdGz.resize(m_lNpt);
    m_vdGx[0] = 0e0; m_vdGy[0] = 0e0; m_vdGz[0] = 0e0;
    for (int64_t lIdxPt = 1; lIdxPt < m_lNpt; ++lIdxPt)
    {
        m_vdGx[lIdxPt] = (m_vdKx[lIdxPt] - m_vdKx[lIdxPt-1])/dDt;
        m_vdGy[lIdxPt] = (m_vdKy[lIdxPt] - m_vdKy[lIdxPt-1])/dDt;
        m_vdGz[lIdxPt] = (m_vdKz[lIdxPt] - m_vdKz[lIdxPt-1])/dDt;
    }
}
