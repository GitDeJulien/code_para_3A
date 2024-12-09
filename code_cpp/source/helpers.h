//
// Created by UTENTE on 13/11/2024.
//

#ifndef SOURCE_HELPERS_H
#define SOURCE_HELPERS_H

#include <string>
#include <iostream>
#include <vector>


double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

#endif //SOURCE_HELPERS_H
