#include "Spiral2D.h"

Spiral2D::Spiral2D() {}

Spiral2D::~Spiral2D() {}

Spiral2D& Spiral2D::operator=(Spiral2D &o)
{
	return *this;
}

void Spiral2D::Update(std::vector<double> *pvdPara)
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

    // k[1]
    dRho = dS*dDtOv*dDtOv;
    dTht = (2e0*M_PI)/(0.5e0)*(dNp/2e0)/(dU)*dRho;
    if (++lCntOv == round(dOv)) {PUSHBACK_KXY(); lCntOv = 0;}

    // k[2] ~ k[inf]
    double dD1Tht = dTht/dDtOv, dD2Tht = 0;
    while (dRho < dKmax)
    {
        dD2Tht = SovD2Tht(dTht, dD1Tht, dNp, dU, dS);
        dD1Tht += dD2Tht*dDtOv;
        dTht += dD1Tht*dDtOv;
        dRho += (0.5e0)/(2e0*M_PI)*(dU)/(dNp/2e0)*dD1Tht*dDtOv;

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

double Spiral2D::SovD2Tht(double dD0Tht, double dD1Tht, double dNp, double dU, double dS)
{
    double dA = 0.25*pow(dU, 2)*(pow(dD0Tht, 2) + 1)/(pow(M_PI, 2)*pow(dNp, 2));
    double dB = 0.5*dD0Tht*pow(dD1Tht, 2)*pow(dU, 2)/(pow(M_PI, 2)*pow(dNp, 2));
    double dC = (0.25*pow(dD0Tht, 2)*pow(dD1Tht, 4)*pow(dU, 2) + 1.0*pow(dD1Tht, 4)*pow(dU, 2) - 1.0*pow(M_PI, 2)*pow(dNp, 2)*pow(dS, 2))/(pow(M_PI, 2)*pow(dNp, 2));

    return SovQuadEq(dA, dB, dC);
}
