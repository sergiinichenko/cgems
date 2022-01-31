#ifndef VECTOR_H
#define VECTOR_H

#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#ifdef DEF_TYPE_FLOAT
    typedef float dType;
#else
    typedef double dType;
#endif

#ifdef DEF_TYPE_INT
    typedef unsigned long int iType;
#else
    typedef int iType;
#endif

//std::vector<dType> operator=(const dType &x);
//std::vector<dType> operator=(const bool &x);

//Multiplication
std::vector<std::vector<dType>> operator*(const std::vector<std::vector<dType>>& x, const std::vector<std::vector<dType>>& y);
std::vector<dType> operator*(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<dType> operator*(const std::vector<dType>& x, const dType& y);
std::vector<dType> operator*(const dType& x, const std::vector<dType>& y);

//Devison
std::vector<dType> operator/(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<dType> operator/(const std::vector<dType>& x, const dType& y);
std::vector<dType> operator/(const dType& x, const std::vector<dType>& y);


//Plus
std::vector<dType> operator+(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<dType> operator+(const std::vector<dType>& x, const dType& y);
std::vector<dType> operator+(const dType& x, const std::vector<dType>& y);

//Plus
std::vector<std::string> operator+(const std::vector<std::string>& x, const std::vector<std::string>& y);
std::vector<std::string> operator+(const std::vector<std::string>& x, const std::string& y);
std::vector<std::string> operator+(const std::string& x, const std::vector<std::string>& y);

//Minus
std::vector<dType> operator-(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<dType> operator-(const std::vector<dType>& x, const dType& y);
std::vector<dType> operator-(const dType& x, const std::vector<dType>& y);
std::vector<dType> operator-(const std::vector<dType>& x);

//exp
std::vector<dType> exp(const std::vector<dType>& x);
dType expNew(const dType& x);


//sqrt
std::vector<dType> sqrt(const std::vector<dType>& x);

//log
std::vector<dType> log(const std::vector<dType>& x);

//log10
std::vector<dType> log10(const std::vector<dType>& x);

//abs
std::vector<dType> abs(const std::vector<dType>& x);

//sum
dType sum(const std::vector<unsigned>& x);
std::vector<dType> sum(const std::vector<std::vector<dType>>& x);
dType sum(const std::vector<dType>& x);
dType sum(const std::vector<dType>& x, const dType& val);
dType sum(const std::vector<dType>& X, const std::vector<int>& ind);

//erf / erfc
std::vector<dType> erfc(const std::vector<dType>& x);
std::vector<dType> erf(const std::vector<dType>& x);

// pow
std::vector<dType> pow(const std::vector<dType>& X, const dType& p);
std::vector<dType> pow(const dType& X, const std::vector<dType>& p);

std::vector<bool> operator<(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<bool> operator<(const std::vector<dType>& x, const dType& y);
std::vector<bool> operator<(const dType& x, const std::vector<dType>& y);

std::vector<bool> operator>(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<bool> operator>(const std::vector<dType>& x, const dType& y);
std::vector<bool> operator>(const dType& x, const std::vector<dType>& y);

std::vector<bool> operator==(const std::vector<dType>& x, const std::vector<dType>& y);
std::vector<bool> operator==(const std::vector<dType>& x, const dType& y);
std::vector<bool> operator==(const dType& x, const std::vector<dType>& y);

std::vector<bool> operator==(const dType& x, const std::vector<dType>& y);

void subset(std::vector<dType>& Inp, const std::vector<iType>& ind, std::vector<dType>& Out);
void subset(std::vector<std::string>& Inp, const std::vector<iType>& ind, std::vector<std::string>& Out);
std::vector<dType> subset(std::vector<dType>& X, const std::vector<iType>& ind);
std::vector<std::string> subset(std::vector<std::string>& X, const std::vector<iType>& ind);
std::vector<dType> condsubset(std::vector<dType>& X, const std::vector<bool>& cond);
std::vector<iType> condsubset(std::vector<iType>& X, const std::vector<bool>& cond);

std::vector<bool> operator&&(const std::vector<bool>& x, const std::vector<bool>& y);
void Clear(std::vector<std::vector<dType>>& X);
void Resize(std::vector<std::vector<dType>>& X);

#endif
