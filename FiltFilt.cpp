#include <algorithm>
#include "FiltFilt.h"
#include "EigenUtil.h"

namespace Digital
{
    namespace Core
    {
        void CoreDoubleN(const double* x, uint32_t xSize, const double* a, const double* b,
            const uint32_t order, double* Y, double* Z)
        {
            double Yi, Xi;
            double b0 = b[0];
            double aOrder = a[order];
            double bOrder = b[order];

            uint32_t i = 0;
            uint32_t j;

            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z[0];          // Filtered value
                Y[i] = Yi;
                for (j = 1; j < order; j++) 
                {
                    Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi; // Update conditions
                }
                Z[order - 1] = bOrder * Xi - aOrder * Yi;
                ++i;
            }
        }

        void CoreDouble2(const double* x, uint32_t xSize, const double* a, const double* b,
                        double*  Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double b0 = b[0];
            double b1 = b[1];

            double Z0 = Z[0];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi - a1 * Yi;
                ++i;
            }

            Z[0] = Z0;
        }

        void CoreDouble3(const double* x, uint32_t xSize, const double* a, const double* b,
                        double* Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double a2 = a[2];
            double b0 = b[0];
            double b1 = b[1];
            double b2 = b[2];

            double Z0 = Z[0];
            double Z1 = Z[1];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi + Z1 - a1 * Yi;
                Z1 = b2 * Xi - a2 * Yi;
                ++i;
            }

            Z[0] = Z0;
            Z[1] = Z1;
        }

        void CoreDouble4(const double* x, uint32_t xSize, const double* a, const double* b,
                        double* Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double a2 = a[2];
            double a3 = a[3];
            double b0 = b[0];
            double b1 = b[1];
            double b2 = b[2];
            double b3 = b[3];

            double Z0 = Z[0];
            double Z1 = Z[1];
            double Z2 = Z[2];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi + Z1 - a1 * Yi;
                Z1 = b2 * Xi + Z2 - a2 * Yi;
                Z2 = b3 * Xi - a3 * Yi;
                i++;
            }

            Z[0] = Z0;
            Z[1] = Z1;
            Z[2] = Z2;
        }

        void CoreDouble5(const double* x, uint32_t xSize, const double* a, const double* b,
                        double* Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double a2 = a[2];
            double a3 = a[3];
            double a4 = a[4];
            double b0 = b[0];
            double b1 = b[1];
            double b2 = b[2];
            double b3 = b[3];
            double b4 = b[4];

            double Z0 = Z[0];
            double Z1 = Z[1];
            double Z2 = Z[2];
            double Z3 = Z[3];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi + Z1 - a1 * Yi;
                Z1 = b2 * Xi + Z2 - a2 * Yi;
                Z2 = b3 * Xi + Z3 - a3 * Yi;
                Z3 = b4 * Xi - a4 * Yi;
                ++i;
            }

            Z[0] = Z0;
            Z[1] = Z1;
            Z[2] = Z2;
            Z[3] = Z3;
        }

        void CoreDouble6(const double* x, uint32_t xSize, const double* a, const double* b,
                        double* Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double a2 = a[2];
            double a3 = a[3];
            double a4 = a[4];
            double a5 = a[5];
            double b0 = b[0];
            double b1 = b[1];
            double b2 = b[2];
            double b3 = b[3];
            double b4 = b[4];
            double b5 = b[5];

            double Z0 = Z[0];
            double Z1 = Z[1];
            double Z2 = Z[2];
            double Z3 = Z[3];
            double Z4 = Z[4];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi + Z1 - a1 * Yi;
                Z1 = b2 * Xi + Z2 - a2 * Yi;
                Z2 = b3 * Xi + Z3 - a3 * Yi;
                Z3 = b4 * Xi + Z4 - a4 * Yi;
                Z4 = b5 * Xi - a5 * Yi;
                ++i;
            }

            Z[0] = Z0;
            Z[1] = Z1;
            Z[2] = Z2;
            Z[3] = Z3;
            Z[4] = Z4;
        }

        void CoreDouble7(const double* x, uint32_t xSize, const double* a, const double* b,
                        double* Y, double* Z)
        {
            double Xi, Yi;
            double a1 = a[1];
            double a2 = a[2];
            double a3 = a[3];
            double a4 = a[4];
            double a5 = a[5];
            double a6 = a[6];
            double b0 = b[0];
            double b1 = b[1];
            double b2 = b[2];
            double b3 = b[3];
            double b4 = b[4];
            double b5 = b[5];
            double b6 = b[6];

            double Z0 = Z[0];
            double Z1 = Z[1];
            double Z2 = Z[2];
            double Z3 = Z[3];
            double Z4 = Z[4];
            double Z5 = Z[5];

            uint32_t i = 0;
            while(i < xSize)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z0;
                Y[i] = Yi;
                Z0 = b1 * Xi + Z1 - a1 * Yi;
                Z1 = b2 * Xi + Z2 - a2 * Yi;
                Z2 = b3 * Xi + Z3 - a3 * Yi;
                Z3 = b4 * Xi + Z4 - a4 * Yi;
                Z4 = b5 * Xi + Z5 - a5 * Yi;
                Z5 = b6 * Xi - a6 * Yi;
                ++i;
            }

            Z[0] = Z0;
            Z[1] = Z1;
            Z[2] = Z2;
            Z[3] = Z3;
            Z[4] = Z4;
            Z[5] = Z5;
        }

        void CoreDoubleNR(const double* x, uint32_t xSize, const double* a, const double* b,
                        uint32_t order, double* Y, double* Z)
        {
            double Yi, Xi;
            double b0 = b[0];
            double aOrder = a[order];
            double bOrder = b[order];

            uint32_t i = xSize - 1;    
            uint32_t j;    
            while(true)
            {
                Xi = x[i];
                Yi = b0 * Xi + Z[0];
                Y[i] = Yi;
                for (j = 1; j < order; j++)
                {
                    Z[j - 1] = b[j] * Xi + Z[j] - a[j] * Yi;
                }
                Z[order - 1] = bOrder * Xi - aOrder * Yi;
                if(i == 0)
                {
                    break;
                }
                --i;
            }
        }
    };

    bool Filter::operator()(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x,
                    const std::vector<double>& z, std::vector<double>& Y, std::vector<double>& Z, bool reverse)
    {
        Y.resize(x.size());
        Z = z;

        if(a.size() != b.size() || a[0] != 1)
        {
            return false;
        }

        size_t order = a.size()-1;
        if(z.size() != order)
        {
            return false;
        }

        if (reverse)
        {
            Core::CoreDoubleNR(x.data(), x.size(), a.data(), b.data(), order, Y.data(), Z.data());
        }
        else
        {
            switch (order) {
                case 1:   Core::CoreDouble2(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                case 2:   Core::CoreDouble3(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                case 3:   Core::CoreDouble4(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                case 4:   Core::CoreDouble5(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                case 5:   Core::CoreDouble6(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                case 6:   Core::CoreDouble7(x.data(), x.size(), a.data(), b.data(), Y.data(), Z.data());  break;
                default:  Core::CoreDoubleN(x.data(), x.size(), a.data(), b.data(), order, Y.data(), Z.data());
            }
        }

        return true;
    }

    FiltFilt::EResult FiltFilt::operator()(std::vector<double>& b, std::vector<double>& a, const std::vector<double>& x, std::vector<double>& Y)
    {
        Filter filter;
        Y.clear();

        if(b.size() > a.size())
            a.resize(b.size(), 0);
        else if(a.size() > b.size())
            b.resize(a.size(), 0);

        if (a[0] == 0.0)
        {
            return EResult::StartsWithZero;
        }

        if (a[0] != 1.0)
        {       
            std::transform(a.begin(), a.end(), a.begin(), [&a](double v) { return v / a[0]; });
            std::transform(b.begin(), b.end(), b.begin(), [&a](double v) { return v / a[0]; });
        }

        size_t Order = a.size();
        size_t nEdge = 3 * (Order - 1);
        if(x.size() <= nEdge)
        {
            return EResult::InputBelowEdge;;
        }

        size_t xLen = x.size();

        size_t kSize = Order-1;
        Eigen::MatrixXd eK = Eigen::MatrixXd::Identity(kSize, kSize);
        for(int i = 0; i < kSize; ++i)
        {
            eK(i, 0) = a[i+1];
            if(i > 0)
            {
                eK(i-1, i) = -1;
            }
        }
        eK(0, 0) += 1;

        Eigen::VectorXd eICTmp1 = Util::toEigen(b, 1);
        Eigen::VectorXd eICTmp2 = Util::toEigen(a, 1);

        eICTmp2 *= b[0];
        eICTmp1 -= eICTmp2;

        Eigen::HouseholderQR<Eigen::MatrixXd> qr(eK);
        Eigen::VectorXd eIC = qr.solve(eICTmp1);
        //slightly slower eIC = eK.colPivHouseholderQr().solve(eICTmp1);

        double x1_2 = 2 * x[0];
        double xf_2 = 2 * x[xLen-1];
    
        std::vector<double> Xi(nEdge, x1_2);
        std::vector<double> Xf(nEdge, xf_2);
        for(int i = 0; i < nEdge; ++i)
        {
            Xi[i] -= x[nEdge-i];
            Xf[i] -= x[(xLen-2)-i];
        }

        std::vector<double> ICxXi1(eIC.begin(), eIC.end());
        double ICxXi1Multiplier = Xi[0];
        std::transform(ICxXi1.begin(), ICxXi1.end(), ICxXi1.begin(),
                [ICxXi1Multiplier](double value) -> double { return value * ICxXi1Multiplier; });

        // #Filter initial reflected signal
        std::vector<double> dum;
        std::vector<double> Zi;
        filter(b, a, Xi, ICxXi1, dum, Zi);

        // #Use the final conditions of the initial part for the actual signal
        std::vector<double> Ys;
        std::vector<double> Zs;
        // #"s"teady state
        filter(b, a, x, Zi, Ys, Zs);

        // #"f"inal conditions
        std::vector<double> Yf;
        filter(b, a, Xf, Zs, Yf, dum);

        // #Filter signal again in reverse order
        std::vector<double> ICxYFnEdge(eIC.begin(), eIC.end());
        double ICxYFnEdgeMultiplier = Yf[nEdge-1];
        std::transform(ICxYFnEdge.begin(), ICxYFnEdge.end(), ICxYFnEdge.begin(),
                [ICxYFnEdgeMultiplier](double value) -> double { return value * ICxYFnEdgeMultiplier; });

        std::vector<double> Zf;
        filter(b, a, Yf, ICxYFnEdge, dum, Zf, true);

        filter(b, a, Ys, Zf, Y, dum, true);
    }
}