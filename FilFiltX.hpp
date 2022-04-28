#include <cstdint>
#include <vector>
#include <algorithm>
#include <iostream>
#include "FilterX.hpp"
#include "eigen-3.4.0/Eigen/Dense"

//#define DEBUG_MODE 1

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}

Eigen::VectorXd toEigen(const vectord& other)
{
    Eigen::VectorXd vec(other.size());
    for(int i = 0; i < other.size(); i++)
        vec(i) = other[i];

    return vec;
}

Eigen::MatrixXd toEigen(const vectord& other, uint32_t columns, uint32_t rows)
{
    Eigen::MatrixXd mx(rows, columns);
    if(rows*columns == other.size())
    {
        for(auto i = 0; i < columns; i++)
        {
            for(int j = 0; j < rows; j++)
            {
                mx(j, i) = other[j+(rows*i)];
            }
        }
    }

    return mx;
}    

vectord fromEigen(const Eigen::VectorXd& other)
{
    vectord vec(other.begin(), other.end());
    return vec;
}

void FiltFiltX(vectord& b, vectord& a, vectord& X, vectord& outY)
{
    outY.clear();

    if(b.size() > a.size())
        a.resize(b.size(), 0);
    else if(a.size() > b.size())
        b.resize(a.size(), 0);

    if (a[0] != 1.0)
    {       
        std::transform(a.begin(), a.end(), a.begin(), [&a](double v) { return v / a[0]; });
        std::transform(b.begin(), b.end(), b.begin(), [&a](double v) { return v / a[0]; });
    }

    size_t Order = a.size();
    size_t nEdge = 3 * (Order - 1);
    if(X.size() <= nEdge)
    {
        std::cout << "X.size() <= nEdge" << std::endl;
        return;
    }

    size_t xLen = X.size();

#ifdef DEBUG_MODE
    std::cout 
    << "X:\n" << X << "\n"
    << "a:\n" << a << "\n"
    << "b:\n" << b << "\n"
    << "Order:\n" << Order << "\n"
    << "nEdge:\n" << nEdge << "\n"
    << "xLen:\n" << xLen << "\n"
    << std::endl;
#endif

    size_t kSize = Order-1;
    vectord K(kSize*kSize, 0.0);

    size_t column = 0;
    while(column < kSize)
    {
        if(column == 0)
        {
            int i = 1;
            while(i < a.size())
            {
                K[i-1] = a[i];
                ++i;
            }

            ++K[0];
        }
        else
        {
            K[(column*kSize)+column] = 1.0;
            K[(column*kSize)+column-1] = -1.0;
        }
        ++column;
    }

#ifdef DEBUG_MODE
    std::cout 
    << "K:\n" << K << "\n"
    << std::endl;
#endif
    

    Eigen::MatrixXd eK = toEigen(K, kSize, kSize);

#ifdef DEBUG_MODE
    std::cout 
    << "eK:\n" << eK << "\n"
    << std::endl;
#endif


	Eigen::VectorXd eICTmp1 = toEigen(vectord(std::next(b.begin()), b.end()));
    Eigen::VectorXd eICTmp2 = toEigen(vectord(std::next(a.begin()), a.end()));

    eICTmp2 *= b[0];
    eICTmp1 -= eICTmp2;

    Eigen::VectorXd eIC(eICTmp1.size());

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(eK);
	eIC = qr.solve(eICTmp1);
    //slightly slower eIC = eK.colPivHouseholderQr().solve(eICTmp1);

#ifdef DEBUG_MODE
    std::cout 
    << "eIC:\n" << eIC << "\n"
    << std::endl;
#endif

    double x1_2 = 2 * X[0];
    double xf_2 = 2 * X[xLen-1];
 
    vectord Xi(nEdge, x1_2);
    for(int i = 0; i < nEdge; i++)
    {
        Xi[i] -= X[nEdge-i];
    }

    vectord Xf(nEdge, xf_2);
    for(int i = 0; i < nEdge; i++)
    {
        Xf[i] -= X[(xLen-2)-i];
    }

#ifdef DEBUG_MODE
    std::cout 
    << "Xi:\n" << Xi << "\n"
    << "Xf:\n" << Xf << "\n"
    << std::endl;
#endif

    //  #1
    vectord ICxXi1 = fromEigen(eIC);
    std::transform(ICxXi1.begin(), ICxXi1.end(), ICxXi1.begin(), [&Xi](double value)
    {
        return value * Xi[0];
    });

#ifdef DEBUG_MODE
    std::cout 
    << "IC * Xi(1, :):\n" << ICxXi1 << "\n"
    << std::endl;
#endif

    vectord dum;
    vectord Zi;
    FilterX(b, a, Xi, ICxXi1, false, dum, Zi);

#ifdef DEBUG_MODE
    std::cout 
    << "dum:\n" << dum << "\n"
    << "Zi:\n" << Zi << "\n"
    << std::endl;
#endif

    // #2
    vectord Ys;
    vectord Zs;
    FilterX(b, a, X, Zi, false, Ys, Zs);

#ifdef DEBUG_MODE
    std::cout 
    << "Ys:\n" << Ys << "\n"
    << "Zs:\n" << Zs << "\n"
    << std::endl;
#endif

    // #3
    vectord Yf;
    FilterX(b, a, Xf, Zs, false, Yf, dum);

#ifdef DEBUG_MODE
    std::cout 
    << "Yf:\n" << Yf << "\n"
    << std::endl;
#endif

    // #4
    vectord ICxYFnEdge = fromEigen(eIC);
    std::transform(ICxYFnEdge.begin(), ICxYFnEdge.end(), ICxYFnEdge.begin(), [nEdge, &Yf](double value)
    {
        return value * Yf[nEdge-1];
    });

#ifdef DEBUG_MODE
    std::cout 
    << "IC * Yf(nEdge, :):\n" << ICxYFnEdge << "\n"
    << std::endl;
#endif

    vectord Zf;
    FilterX(b, a, Yf, ICxYFnEdge, true, dum, Zf);

#ifdef DEBUG_MODE
    std::cout 
    << "Zf:\n" << Zf << "\n"
    << std::endl;
#endif

    // #5 final
    FilterX(b, a, Ys, Zf, true, outY, dum);

#ifdef DEBUG_MODE
    std::cout 
    << "Y:\n" << outY << "\n"
    << std::endl;
#endif
}