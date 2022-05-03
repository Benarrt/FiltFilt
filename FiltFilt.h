#ifndef FILT_FITL_H
#define FILT_FITL_H

#include <vector>

namespace Digital
{
    class FiltFilt
    {
        public:
        enum class EResult : unsigned
        {
            Succes = 0,
            StartsWithZero,
            InputBelowEdge,
            UnkownError
        };

        /*! \brief Performs zero-phase digital filtering. Filtering is done by first procesing data forward, then backward.
            \param b - The numerator coefficient vector of the filter. Will be padded with 0's to the lenght of a if shorter than a.
            \param a - The denominator coefficient vector of the filter. Will be added with 0's to the lenght of b if shorter than b. 
            If a[0] is not 1 then a and b get normalized by a[0], and a[0] is set to 1. a[0] can not be 0.
            \param x - Input data to be filtered.
            \param Y - Filtered output data.
            \return EFiltFiltResult::Succes - Filtering was finished with success,
            EFiltFiltResult::StartsWithZero - Error, filtering not finished A[0] is equal 0,
            EFiltFiltResult::InputBelowEdge - Error, x is shorter that the the edge ((max(b, a)-1)*3),
            EFiltFiltResult::UnkownError - Error during filtering process, output Y may contain some invalid data.
        */
        EResult operator()(std::vector<double>& b, std::vector<double>& a, const std::vector<double>& x, std::vector<double>& Y);
    };

    class Filter
    {
        public:
        /*! \brief Filters data.
            \param b - The numerator coefficient vector of the filter. Must be the size of a.
            \param a - The denominator coefficient vector of the filter. Must be the size of b. First element must equal to 1.
            \param x - Input data to be filtered.
            \param z - Initial conditions, size must be a.size-1 .
            \param Y - Filter output data.
            \param Z - Output final conditions.
            \param reverse - True if filter should process in backward direction.
            \return true - Filtering was finished with success,
            false - Error during filtering process, output Y may contain some invalid data.
        */
        bool operator()(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x,
                        const std::vector<double>& z, std::vector<double>& Y, std::vector<double>& Z, bool reverse = false);
    };
    
} // namespace Digital


#endif //FILT_FITL_H