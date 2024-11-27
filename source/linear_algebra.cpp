#ifndef _LINEAR_ALGEBRA_CPP
#define _LINEAR_ALGEBRA_CPP

#include "linear_algebra.h"


LinearAlgebra::LinearAlgebra(){};

std::vector<double> LinearAlgebra::LU(const Matrix& A, const std::vector<double> b){

    int N = b.size();

    Matrix L(N, N);
    Matrix U(N, N);

    U = A;

    std::vector<double> X;
    X.resize(N);

    if (A.getCols() != N || A.getCols() != A.getCols()) {
        std::cerr << "Error: dimensions don't agree in LU function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        L(i,i) = 1.;
    }

    //building L and U via gauss elimination
    for (int k = 0; k < N; ++k) {
        for (int i = k+1; i < N; i++) {
            double multiplier = U(i,k)/U(k,k);
            L(i,k) = multiplier;
            U(i,k) = 0.;
            for (int j = k+1; j < N; j++) {
                U(i,j) -= multiplier*U(k,j);
            }
        }
    }

    //backward solving Ux=y
    for (int i = 1; i <= N; ++i) {
        double s = 0.;
        for (int j = 1; j < i; ++j){
            s += U(N-i,N-j)*X[N-j];
        }
        X[N-i] = (b[N-i] - s)/U(N-i,N-i);
    }

    //forward solving Ly = b
    for (int i = 0; i < N; ++i) {
        double s = 0.;
        for (int j = 0; j < i; ++j)
            s += L(i,j)*X[j];
        X[i] = (b[i] - s)/L(i,i);
    }

    return X;

}

Matrix LinearAlgebra::LU(const Matrix& A, const Matrix& B){

    int N = A.getRows();

    Matrix L(N, N);
    Matrix U(N, N);

    U = A;

    Matrix X(A.getRows(),B.getCols());
    Matrix Y(A.getRows(),B.getCols());

    if (A.getCols() != N || A.getCols() != A.getCols()) {
        std::cerr << "Error: dimensions don't agree in LU function!" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        L(i,i) = 1.;
    }

    //building L and U via gauss elimination
    for (int k = 0; k < N; ++k) {
        for (int i = k+1; i < N; i++) {
            double multiplier = U(i,k)/U(k,k);
            L(i,k) = multiplier;
            U(i,k) = 0.;
            for (int j = k+1; j < N; j++) {
                U(i,j) -= multiplier*U(k,j);
            }
        }
    }
    //forward solving Ly = b
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < B.getCols(); ++j) {
            double s = 0.;
            for (int k = 0; k < i; ++k) { // Solo termini precedenti
                s += L(i,k) * Y(k,j);
            }
            Y(i,j) = (B(i,j) - s) / L(i,i);
        }
    }


    //backward solving Ux=y
    for (int i = N-1; i >= 0; --i) {
        for (int j = 0; j < B.getCols(); ++j) {
            double s = 0.;
            for (int k = i+1; k < N; ++k) { // Solo termini successivi
                s += U(i,k) * X(k,j);
            }
            X(i,j) = (Y(i,j) - s) / U(i,i);
        }
    }

    //X.print();
    return X;

}

std::vector<double> LinearAlgebra::BiCGStab(const Matrix& A, const std::vector<double> b, int maxIterations = 1000, double tol = 1e-6){
    
    int N = b.size();

    std::vector<double> X(N, 0.0);             // Initial guess (zero vector)
    std::vector<double> r = b;                 // Residual vector r = b - A * X (initially, b)
    std::vector<double> r0 = r;                // Copy of the initial residual
    std::vector<double> p(N, 0.0), v(N, 0.0), s(N, 0.0), t(N, 0.0), h(N,0.0);
    
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_prev = 0.0, beta = 0.0;
    
    double normB = norm(b);
    if (normB < tol) return X; // If b is small enough, return X = 0
    
    p = r;           
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        rho_prev = rho;
        rho = dot(r0, r);
        beta = (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < N; ++i){
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        if (std::abs(rho_prev) < tol) break;       // Breakdown check
        
        v = A.MatrixVectorProduct(p);
        alpha = rho / dot(r0, v);
        
        for (int i = 0; i < N; ++i){
            h[i] = X[i] + alpha * p[i];
            s[i] = r[i] - alpha * v[i];
        }
        
        if (norm(h) < tol) {
            for (int i = 0; i < N; ++i)
                X[i] += alpha * p[i];
            break;
        }
        
        t = A.MatrixVectorProduct(s);
        omega = dot(t, s) / dot(t, t);
        
        for (int i = 0; i < N; ++i)
            X[i] = h[i] + omega * s[i];
        
        for (int i = 0; i < N; ++i)
            r[i] = s[i] - omega * t[i];
        
        //normB = norm(b);
        if (norm(r) < tol*normB) break;
        
    }

    return X;
}

Matrix LinearAlgebra::BiCGStab(const Matrix& A, const Matrix& B, int maxIterations = 1000, double tol = 1e-6){

    Matrix X(A.getCols(),B.getRows());             // Initial guess (zero vector)
    Matrix r = B;                 // Residual vector r = b - A * X (initially, b)
    Matrix r0 = r;                // Copy of the initial residual
    Matrix p(A.getCols(),B.getRows()), v(A.getCols(),B.getRows()),
    s(A.getCols(),B.getRows()), t(A.getCols(),B.getRows()),
    h(A.getCols(),B.getRows());

    double rho(1.0), alpha(1.0), omega = 1.0;
    double rho_prev = 0.0, beta = 0.0;

    double normB = norm(B);
    if (normB < tol) return X; // If b is small enough, return X = 0

    p = r;

    for (int iter = 0; iter < maxIterations; ++iter) {
        rho_prev = rho;
        rho = norm(r0*r);
        beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p + v*(-omega))*beta;
        if (std::abs(rho_prev) < tol) break;       // Breakdown check

        v = A*p;
        alpha = rho / norm(r0*v);

        h = X + p * alpha;
        s = r + v*(- alpha);

        if (norm(h) < tol) {
            X = X + alpha * p;
            break;
        }

        t = A*s;
        omega = norm(t*s) / norm((t*t));

        X = h + omega * s;

        r = s - omega * t;

        //normB = norm(b);
        if (norm(r) < tol*normB) break;

    }

    return X;
}

Matrix LinearAlgebra::solveLiap(const Matrix & A, const Matrix & B, const Matrix & C, int nit) {
    Matrix X = C;
    Matrix X1 = X;
    Matrix Ia(A.getRows(),A.getCols(),1);
    Matrix Ib(B.getRows(),B.getCols(),1);
    //Coefficients de stabilisation pour la bonne convergence de la méthode
    std::vector<double> coefficients = computeAdiCoefficients(A.getRows(), A(0,0), A(1,0), B(0,0), B(1,0), nit);
    //std::vector<double> coefficients = computeAdiCoefficients(3,2, -1, 2, -1, nit);
    double tol = 1e-7;
    for(int i=0; i<nit; i++){
        X1 = LU(A - coefficients[i]*Ia,X*(B-coefficients[i]*Ib) + C);
        X = trMatrix(LU(trMatrix(B - coefficients[nit-i]*Ib),trMatrix((A-coefficients[nit-i]*Ib)*X1 - C)));
        if(norm(A*X - X*B - C)<tol){
            (A*X - X*B - C).print();
            return X;
        }
    }
    return X;
}



#endif