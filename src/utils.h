// Copyright (2020) Chongliang Luo, Rui Duan, Jiayi Tong and Yong Chen
//     
//     This file is part of pda
//     
//     Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//     
//     http://www.apache.org/licenses/LICENSE-2.0
//     
//     Unless required by applicable law or agreed to in writing, software
//     distributed under the License is distributed on an "AS IS" BASIS,
//     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//     limitations under the License.


// Copyright (C) 2017-2018  Wenjie Wang <wjwang.stat@gmail.com>
//
// This file is part of the R package intsurv.
//
// The R package intsurv is free software: You can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any later
// version (at your option). See the GNU General Public License at
// <https://www.gnu.org/licenses/> for details.
//
// The R package intsurv is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#ifndef UTILS_H
#define UTILS_H

#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

// #include <armadillo>
#include <RcppArmadillo.h>

namespace Intsurv {

    // compare double-precision numbers for almost equality
    inline bool isAlmostEqual(double A, double B)
    {
        double MaxRelDiff {std::numeric_limits<double>::epsilon()};
        // compute the difference.
        double diff = std::abs(A - B);
        A = std::abs(A);
        B = std::abs(B);
        // Find the largest
        double largest = (B > A) ? B : A;
        if (diff <= largest * MaxRelDiff) {
            return true;
        } else {
            return false;
        }
    }

    // function checking if there exists any duplicates
    inline bool any_duplicated(const arma::vec& x)
    {
        std::unordered_set<double> seen;
        bool res {false};
        for (size_t i {0}; i < x.n_rows; ++i) {
            res = ! seen.insert(x(i)).second;
            if (res) break;
        }
        return res;
    }
    // function checking if there exists any duplicates
    inline arma::uvec duplicated(const arma::vec& x, bool fromLast = false)
    {
        std::unordered_set<double> seen;
        std::vector<unsigned int> res;
        if (fromLast) {
            for (size_t i {1}; i <= x.n_rows; ++i) {
                if (! seen.insert(x(x.n_rows - i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        } else {
            for (size_t i {0}; i < x.n_rows; ++i) {
                if (! seen.insert(x(i)).second) {
                    // if duplicated, add index to vector res
                    res.push_back(i);
                }
            }
        }
        return arma::conv_to<arma::uvec>::from(res);
    }
    // function that returns the indices of first unique elements from last
    // inline arma::uvec find_unique_last(const arma::vec& x)
    // {
    //     return std::reverse(x.n_elem - 1 -
    //                          arma::find_unique(std::reverse(x)));
    // }
    
    // set intersection for vector a and vector b
    // armadillo vector has just one template type parameter
    template <typename T, template <typename> class ARMA_VEC_TYPE>
    inline ARMA_VEC_TYPE<T> vec_intersection(const ARMA_VEC_TYPE<T>& a,
                                             const ARMA_VEC_TYPE<T>& b)
    {
        std::vector<T> res;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                              std::back_inserter(res));
        std::reverse(res.begin(), res.end());
        return arma::sort(arma::conv_to<ARMA_VEC_TYPE<T>>::from(res));
    }

    // cumulative sum in possibly reverse order
    inline arma::vec cum_sum(const arma::vec& x,
                             const bool reversely = false)
    {
        // if cumsum reversely
        if (reversely) {
            const unsigned long int n_x {x.n_rows};
            arma::vec res {arma::zeros(n_x)};
            double tmp {0.0};
            for (size_t i {1}; i <= n_x; ++i) {
                tmp += x[n_x - i];
                res[n_x - i] = tmp;
            }
            return res;
        }
        // otherwise, using arma::cumsum
        return arma::cumsum(x);
    }
    // column-wise cumulative sum in possibly reverse order
    inline arma::mat cum_sum_cols(const arma::mat& x,
                                  const bool reversely = false)
    {
        // if cumsum reversely
        if (reversely) {
            const unsigned long int n_x = x.n_rows;
            arma::mat tmp {arma::zeros(1, x.n_cols)};
            arma::mat res {x};
            for (size_t i {1}; i <= n_x; ++i) {
                tmp += x.row(n_x - i);
                res.row(n_x - i) = tmp;
            }
            return res;
        }
        // otherwise, using arma::cumsum
        return arma::cumsum(x, 0);
    }

    // aggregate sum of a vector based on same indices
    inline arma::vec aggregate_sum(const arma::vec& x,
                                   const arma::vec& indices,
                                   const bool simplify = true,
                                   const bool cumulative = false,
                                   const bool reversely = false)
    {
        const unsigned long int n_x {x.size()};
        arma::vec uniInd {arma::unique(indices)};
        const unsigned long int n_uniInd {uniInd.n_rows};
        // the x's having a same index are summed
        arma::vec sumVec {arma::zeros(n_uniInd)};
        for (size_t i {0}; i < n_uniInd; ++i) {
            for (size_t j {0}; j < n_x; ++j) {
                if (isAlmostEqual(uniInd[i], indices[j])) {
                    sumVec[i] += x[j];
                }
            }
        }
        if (cumulative) {
            sumVec = cum_sum(sumVec, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify) {
            return sumVec;
        }
        // else
        arma::vec out {arma::zeros(n_x)};
        for (size_t i {0}; i < n_x; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {     // indices[i] ??
                    out[i] = sumVec[j];
                    break;
                }
            }
        }
        return out;
    }
    // column-wise aggregrate sum
    inline arma::mat aggregate_sum_cols(const arma::mat& x,
                                        const arma::vec& indices,
                                        const bool simplify = true,
                                        const bool cumulative = false,
                                        const bool reversely = false)
    {
        // if it does need aggregate
        const unsigned long int x_nrows {x.n_rows};
        arma::vec uniInd {arma::unique(indices)};
        const unsigned long int n_uniInd {uniInd.n_rows};
        // the x's having a same index are summed
        arma::mat sumMat {arma::zeros(n_uniInd, x.n_cols)};
        for (size_t i {0}; i < n_uniInd; ++i) {
            for (size_t j {0}; j < x_nrows; ++j) {
                if (isAlmostEqual(uniInd[i], indices[j])) {
                    sumMat.row(i) += x.row(j);
                }
            }
        }
        if (cumulative) {
            sumMat = cum_sum_cols(sumMat, reversely);
        }
        // if simplify the sum results to unique and sorted indices
        if (simplify) {
            return sumMat;
        }
        // else
        arma::mat out {arma::zeros(arma::size(x))};
        for (size_t i {0}; i < x_nrows; ++i) {
            for (size_t j {0}; j < n_uniInd; ++j) {
                if (isAlmostEqual(indices[i], uniInd[j])) {
                    out.row(i) = sumMat.row(j);
                    break;
                }
            }
        }
        return out;
    }

    // inline handy functions
    inline arma::vec mat2vec(const arma::mat& x) {
        return arma::conv_to<arma::vec>::from(x);
    }
    inline double vec2num(const arma::vec& x) {
        return arma::as_scalar(x);
    }

    // function template for crossprod of two matrix-like objects
    template <typename T_matrix_like>
    inline arma::mat crossprod(T_matrix_like X, T_matrix_like Y)
    {
        return X.t() * Y;
    }
    template <typename T_matrix_like>
    inline arma::mat crossprod(T_matrix_like X)
    {
        return X.t() * X;
    }
    // function template for tcrossprod of two matrix-like objects
    template <typename T_matrix_like>
    inline arma::mat tcrossprod(T_matrix_like X, T_matrix_like Y)
    {
        return X * Y.t();
    }
    template <typename T_matrix_like>
    inline arma::mat tcrossprod(T_matrix_like X)
    {
        return X * X.t();
    }

    // function that computes L2-norm
    inline double norm(const arma::vec& x)
    {
        return std::sqrt(vec2num(crossprod(x)));
    }
    inline double norm(const arma::vec& x, const arma::vec& y)
    {
        return std::sqrt(vec2num(crossprod(x, y)));
    }

}

#endif
