#include "Spiral2D.h"

Spiral2D::Spiral2D() {}

Spiral2D::~Spiral2D() {}

Spiral2D& Spiral2D::operator=(Spiral2D &o)
{
	return *this;
}

void Spiral2D::Update(double dKtht1, double dKtht2, double dTht0, double dKmax, double dS, double dDt, double dOv)
{
    std::list<double> ldKx(0), ldKy(0);
    double dDtOv = dDt/dOv;
    int64_t lCntOv = 0;
#define PUSHBACK_KXY() \
    {\
        ldKx.push_back(dRho*cos(dTht + dTht0));\
        ldKy.push_back(dRho*sin(dTht + dTht0));\
    }

    // k[0]
    double dRho = 0;
    double dTht = 0;
    PUSHBACK_KXY(); lCntOv = 0;

    // k[1] ~ k[inf]
    double dD1Tht = 0, dD2Tht = 0;
    while (dRho < dKmax)
    {
        dD2Tht = SovD2Tht(dTht, dD1Tht, dKtht1, dKtht2, dS);
        dD1Tht += dD2Tht*dDtOv;
        dTht += dD1Tht*dDtOv;
        dRho = dKtht1*dTht + dKtht2*dTht*dTht;

        if (++lCntOv == round(dOv)) {PUSHBACK_KXY(); lCntOv = 0;}
    }
#undef PUSHBACK_KXY
    m_lNpt = ldKx.size();
    
    // copy to member variables
    m_vdKx.assign(ldKx.begin(), ldKx.end());
    m_vdKy.assign(ldKy.begin(), ldKy.end());
    m_vdKz.clear();

    // update gradient vector
    m_vdGx.resize(m_lNpt); m_vdGy.resize(m_lNpt); m_vdGz.clear();
    m_vdGx[0] = 0e0; m_vdGy[0] = 0e0;
    for (int64_t lIdxPt = 1; lIdxPt < m_lNpt; ++lIdxPt)
    {
        m_vdGx[lIdxPt] = (m_vdKx[lIdxPt] - m_vdKx[lIdxPt-1])/dDt;
        m_vdGy[lIdxPt] = (m_vdKy[lIdxPt] - m_vdKy[lIdxPt-1])/dDt;
    }
}

double Spiral2D::SovD2Tht(double dD0Tht, double dD1Tht, double dKtht1, double dKtht2, double dS)
{
    // double dA = 0.25*pow(dU, 2)*(pow(dD0Tht, 2) + 1)/(pow(M_PI, 2)*pow(dNp, 2));
    // double dB = 0.5*dD0Tht*pow(dD1Tht, 2)*pow(dU, 2)/(pow(M_PI, 2)*pow(dNp, 2));
    // double dC = (0.25*pow(dD0Tht, 2)*pow(dD1Tht, 4)*pow(dU, 2) + 1.0*pow(dD1Tht, 4)*pow(dU, 2) - 1.0*pow(M_PI, 2)*pow(dNp, 2)*pow(dS, 2))/(pow(M_PI, 2)*pow(dNp, 2));
    double dA = pow(dD0Tht, 4)*pow(dKtht2, 2) + 2*pow(dD0Tht, 3)*dKtht1*dKtht2 + pow(dD0Tht, 2)*pow(dKtht1, 2) + 4*pow(dD0Tht, 2)*pow(dKtht2, 2) + 4*dD0Tht*dKtht1*dKtht2 + pow(dKtht1, 2);
    double dB = 2*pow(dD1Tht, 2)*(2*pow(dD0Tht, 3)*pow(dKtht2, 2) + 3*pow(dD0Tht, 2)*dKtht1*dKtht2 + dD0Tht*pow(dKtht1, 2) + 4*dD0Tht*pow(dKtht2, 2) + 2*dKtht1*dKtht2);
    double dC = pow(dD0Tht, 4)*pow(dD1Tht, 4)*pow(dKtht2, 2) + 2*pow(dD0Tht, 3)*pow(dD1Tht, 4)*dKtht1*dKtht2 + pow(dD0Tht, 2)*pow(dD1Tht, 4)*pow(dKtht1, 2) + 12*pow(dD0Tht, 2)*pow(dD1Tht, 4)*pow(dKtht2, 2) + 12*dD0Tht*pow(dD1Tht, 4)*dKtht1*dKtht2 + 4*pow(dD1Tht, 4)*pow(dKtht1, 2) + 4*pow(dD1Tht, 4)*pow(dKtht2, 2) - pow(dS, 2);

    return SovQuadEq(dA, dB, dC);
}
