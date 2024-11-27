#ifndef _MATRIX_H
#define _MATRIX_H

#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "data.h"

class Matrix {

private:

    std::vector<double> data;
    int rows, cols;

public:

    // Constructor prototype
    /*
    Matrix(int r,int c): rows(r), cols(c){
        std::vector<double> data;
        for(size_t i=0;i<rows;i++){
            for(size_t j=0;j<cols;j++){
                data[i*cols+j] = 0;
            }
        }
    };
    */
    Matrix(int r, int c, int type) : rows(r), cols(c), data(r*c, 0){
        for(size_t i=0;i<rows;i++){
            for(size_t j=0;j<cols;j++){
                if(i == j && type == 1){
                    data[i+cols*j] = 1;
                }
            }
        }
    };

    Matrix(int r, int c) : rows(r), cols(c), data(r * c, 0){};

    double& operator()(int i, int j) { return data[i * cols + j]; };
    const double& operator()(int i, int j) const { return data[i * cols + j]; };
    std::vector<double> operator()(int i) const { //Extraire une ligne
        std::vector<double> row;
        for(size_t j=0; j<cols; j++){
            row.push_back(data[i*cols + j]);
        }
        return row;
    };

    // Method to multiply a matrix with a vector
    std::vector<double> MatrixVectorProduct(const std::vector<double>& x) const;

    // Method to multiply two matrices
    Matrix MatrixMatrixProduct(const Matrix& B) const;

    // Method to add two matrices
    Matrix AddMatrix(const Matrix& B) const;

    // Method to multiply matrix by a scalar value
    Matrix ScalarMultiply(double lambda) const;

    // Static method to generate an identity matrix
    static Matrix Identity(int size);

    void print() const;

    int getRows() const;

    int getCols() const;
};


#endif