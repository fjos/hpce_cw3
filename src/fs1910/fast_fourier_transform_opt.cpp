#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>
#include <thread>
#include "tbb/parallel_for.h"
#include "tbb/task_group.h"

namespace hpce
{
namespace fs1910
{

class fast_fourier_transform_opt
    : public fourier_transform
{
private:
    size_t getNumThreads()
    {
        // size_t x = std::thread::hardware_concurrency();
        // if (x == 0)
        // {
        //     x = 1 << 22;
        // }else if (x ==1)
        // {
        //     x = 1 << 22
        // } else if (x == 2)
        // {
        //     x = 
        // }

        return 1 << 1;
    }   
    const size_t K = getNumThreads();
    const size_t group_size = 1 << 1;
protected:
    /* Standard radix-2 FFT only supports binary power lengths */
    virtual size_t calc_padded_size(size_t n) const
    {
        assert(n != 0);

        size_t ret = 1;
        while (ret < n)
        {
            ret <<= 1;
        }

        return ret;
    }

    virtual void forwards_impl(
        size_t n,   const std::complex<double> &wn,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const
    {
        assert(n > 0);

        if (n == 1)
        {
            pOut[0] = pIn[0];
        }
        else if (n == 2)
        {
            pOut[0] = pIn[0] + pIn[sIn];
            pOut[sOut] = pIn[0] - pIn[sIn];
        }
        else if (n == 4)
        {
            forwards_impl(2, wn * wn, pIn, 2 * sIn, pOut, sOut);
            forwards_impl(2, wn * wn, pIn + sIn, 2 * sIn, pOut + sOut * 2, sOut);
            // /for j = 1, 2
            std::complex<double> w = std::complex<double>(1.0, 0.0);

            std::complex<double> t1 = w * pOut[2];
            std::complex<double> t2 = pOut[0] - t1;
            pOut[0] = pOut[0] + t1;
            pOut[2] = t2;
            w = w * wn;

            t1 = w * pOut[3];
            t2 = pOut[1] - t1;
            pOut[1] = pOut[1] + t1;
            pOut[3] = t2;
            // w = w * wn;

        }
        else
        {

            size_t m = n / 2;

            if (n > group_size)
            {
                tbb::task_group group;

                group.run([&]()
                {
                    forwards_impl(m, wn * wn, pIn, 2 * sIn, pOut, sOut);
                });
                group.run([&]()
                {
                    forwards_impl(m, wn * wn, pIn + sIn, 2 * sIn, pOut + sOut * m, sOut);
                });

                group.wait();
            }
            else
            {
                forwards_impl(m, wn * wn, pIn, 2 * sIn, pOut, sOut);
                forwards_impl(m, wn * wn, pIn + sIn, 2 * sIn, pOut + sOut * m, sOut);
            }

            if (m > K)
            {
                tbb::parallel_for<size_t>( (size_t)0, (m / K), [ = ](size_t j0)
                {
                    std::complex<double> w = std::complex<double>(1.0, 0.0);
                    w = std::pow(wn, j0 * K);
                    for (size_t j1 = 0; j1 < K; j1++)
                    {
                        size_t j = j0 * K + j1;
                        std::complex<double> t1 = w * pOut[m + j];
                        std::complex<double> t2 = pOut[j] - t1;
                        pOut[j] = pOut[j] + t1;               /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
                        pOut[j + m] = t2;                        /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
                        w = w * wn;
                    }
                });
            }
            else
            {
                std::complex<double> w = std::complex<double>(1.0, 0.0);
                for (size_t j = 0; j < m; j++)
                {
                    std::complex<double> t1 = w * pOut[m + j];
                    std::complex<double> t2 = pOut[j] - t1;
                    pOut[j] = pOut[j] + t1;               /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
                    pOut[j + m] = t2;                        /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
                    w = w * wn;
                }
            }

        }
    }

    virtual void backwards_impl(
        size_t n,   const std::complex<double> &wn,
        const std::complex<double> *pIn, size_t sIn,
        std::complex<double> *pOut, size_t sOut
    ) const
    {
        complex_t reverse_wn = 1.0 / wn;
        forwards_impl(n, reverse_wn, pIn, sIn, pOut, sOut);

        double scale = 1.0 / n;
        for (size_t i = 0; i < n; i++)
        {
            pOut[i] = pOut[i] * scale;
        }
    }

public:
    virtual std::string name() const
    {
        return "hpce.fs1910.fast_fourier_transform_opt";
    }

    virtual bool is_quadratic() const
    {
        return false;
    }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_opt()
{
    return std::make_shared<fast_fourier_transform_opt>();
}
}; // namespace fs1910
}; // namespace hpce
