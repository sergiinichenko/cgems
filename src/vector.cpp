#include "vector.h"

// Multiply
std::vector<std::vector<dType>> operator*(const std::vector<std::vector<dType>>& x, const std::vector<std::vector<dType>>& y) {
    std::vector<std::vector<dType>> res(x.size());
    for (int i = 0; i < x.size(); i++) {
        res[i].resize(x[i].size());
        std::transform(x[i].begin(),        x[i].end(),        y[i].begin() ,  res[i].begin(),     std::multiplies<dType>());
    }
    return res;
}

std::vector<dType> operator*(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::multiplies<dType>());
    return res;
}

std::vector<dType> operator*(const std::vector<dType>& x, const dType& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::multiplies<dType>(), y));
    return res;
}

std::vector<dType> operator*(const dType& x, const std::vector<dType>& y) {
    std::vector<dType> res(y.size());
    std::transform(y.begin(),       y.end(),       res.begin(),     std::bind2nd(std::multiplies<dType>(), x));
    return res;
}


// Divide
std::vector<dType> operator/(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::divides<dType>());
    return res;
}

std::vector<dType> operator/(const std::vector<dType>& x, const dType& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::divides<dType>(), y));
    return res;
}

std::vector<dType> operator/(const dType& x, const std::vector<dType>& y) {
    std::vector<dType> res(y.size());
    std::vector<dType> xx(y.size(), x);
    std::transform(xx.begin(),        xx.end(),        y.begin() ,  res.begin(),     std::divides<dType>());
    return res;
}


// Plus
std::vector<dType> operator+(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),         x.end(),         y.begin(),     res.begin(),    std::plus<dType>());
    return res;
}

std::vector<dType> operator+(const std::vector<dType>& x, const dType& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),         x.end(),         res.begin(),    std::bind2nd(std::plus<dType>(), y));
    return res;
}

std::vector<dType> operator+(const dType& x, const std::vector<dType>& y) {
    std::vector<dType> res(y.size());
    std::transform(y.begin(),         y.end(),         res.begin(),    std::bind2nd(std::plus<dType>(), x));
    return res;
}


// Plus
std::vector<std::string> operator+(const std::vector<std::string>& x, const std::vector<std::string>& y) {
    std::vector<std::string> res(x.size());
    std::transform(x.begin(),         x.end(),         y.begin(),     res.begin(),    std::plus<std::string>());
    return res;
}

std::vector<std::string> operator+(const std::vector<std::string>& x, const std::string& y) {
    std::vector<std::string> res(x.size());
    std::transform(x.begin(),         x.end(),         res.begin(),    std::bind2nd(std::plus<std::string>(), y));
    return res;
}

std::vector<std::string> operator+(const std::string& x, const std::vector<std::string>& y) {
    std::vector<std::string> res(y.size());
    std::transform(y.begin(),         y.end(),         res.begin(),    std::bind2nd(std::plus<std::string>(), x));
    return res;
}


// Minus
std::vector<dType> operator-(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),         x.end(),         y.begin(),     res.begin(),    std::minus<dType>());
    return res;
}

std::vector<dType> operator-(const std::vector<dType>& x, const dType& y) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),         x.end(),         res.begin(),    std::bind2nd(std::minus<dType>(), y));
    return res;
}

std::vector<dType> operator-(const dType& x, const std::vector<dType>& y) {
    std::vector<dType> res(y.size());
    std::transform(y.begin(),       y.end(),         res.begin(),    std::bind2nd(std::minus<dType>(), x));
    std::transform(res.begin(),     res.end(),       res.begin(),    std::bind2nd(std::multiplies<dType>(), -1.0));
    return res;
}

std::vector<dType> operator-(const std::vector<dType>& x) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::multiplies<dType>(), -1.0));
    return res;
}

std::vector<dType> exp(const std::vector<dType>& x) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     static_cast<dType (*)(dType)>(std::exp));
    return res;
}

dType expNew(const dType& x) {
    int d;
    dType f, e = 2.7182818285, res = 1.0, re = 0.3678794412;
    d = std::floor(x);
    f = x - d;
    if (x > 0) {
        for (int i = 0; i < d; i++) {
            res *= e;
        }
    } else {
        d *= -1.0;
        for (int i = 0; i < d; i++) {
            res *= re;
        }
    }
    f =1.0 + f*(1.0 + f*0.5*(1.0 + f * 0.33333333333*(1.0 + f*0.25*(1.0 + f * 0.2))));
    res *= f;
    return res;
}

std::vector<dType> abs(const std::vector<dType>& x) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     static_cast<dType (*)(dType)>(std::fabs));
    return res;
}

std::vector<dType> sqrt(const std::vector<dType>& x) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     static_cast<dType (*)(dType)>(std::sqrt));
    return res;
}

std::vector<dType> log(const std::vector<dType>& x) {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     static_cast<dType (*)(dType)>(std::log));
    return res;
}

std::vector<dType> log10(const std::vector<dType>& x)  {
    std::vector<dType> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     static_cast<dType (*)(dType)>(std::log10));
    return res;
}




// Less
std::vector<bool> operator<(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::less<dType>());
    return res;
}

std::vector<bool> operator<(const std::vector<dType>& x, const dType& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::less<dType>(), y));
    return res;
}

std::vector<bool> operator<(const dType& x, const std::vector<dType>& y) {
    std::vector<bool> res(y.size());
    std::transform(y.begin(),       y.end(),       res.begin(),     std::bind2nd(std::less<dType>(), x));
    return res;
}



// greater
std::vector<bool> operator>(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::greater<dType>());
    return res;
}

std::vector<bool> operator>(const std::vector<dType>& x, const dType& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::greater<dType>(), y));
    return res;
}

std::vector<bool> operator>(const dType& x, const std::vector<dType>& y) {
    std::vector<bool> res(y.size());
    std::transform(y.begin(),       y.end(),       res.begin(),     std::bind2nd(std::greater<dType>(), x));
    return res;
}



// Equal
std::vector<bool> operator==(const std::vector<dType>& x, const std::vector<dType>& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::equal_to<dType>());
    return res;
}

std::vector<bool> operator==(const std::vector<dType>& x, const dType& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),       x.end(),       res.begin(),     std::bind2nd(std::equal_to<dType>(), y));
    return res;
}

std::vector<bool> operator==(const dType& x, const std::vector<dType>& y) {
    std::vector<bool> res(y.size());
    std::transform(y.begin(),       y.end(),       res.begin(),     std::bind2nd(std::equal_to<dType>(), x));
    return res;
}





std::vector<bool> operator&&(const std::vector<bool>& x, const std::vector<bool>& y) {
    std::vector<bool> res(x.size());
    std::transform(x.begin(),        x.end(),        y.begin() ,  res.begin(),     std::logical_and<bool>());
    return res;
}



void subset(std::vector<dType>& Inp, const std::vector<iType>& ind, std::vector<dType>& Out) {
    int j;
    Out.reserve(ind.size());
	for (int i = 0; i < ind.size(); i++){
        j = ind[i];
        Out.push_back( Inp[j] );
    }
}

void subset(std::vector<std::string>& Inp, const std::vector<iType>& ind, std::vector<std::string>& Out) {
    int j;
    Out.reserve(ind.size());
    for (int i = 0; i < ind.size(); i++){
        j = ind[i];
        Out.push_back( Inp[j] );
    }
}

std::vector<dType> subset(std::vector<dType>& X, const std::vector<int>& ind) {
    int j;
    //std::vector<dType> res(ind.size());
    std::vector<dType> res;
    res.reserve(ind.size());
    for (int i = 0; i < ind.size(); i++){
        j = ind[i];
        res.push_back( X[j] );
    }
    return res;
}

std::vector<std::string> subset(std::vector<std::string>& X, const std::vector<int>& ind) {
    int j;
    //std::vector<std::string> res(ind.size());
    std::vector<std::string> res;
    res.reserve(ind.size());
    for (int i = 0; i < ind.size(); i++){
        j = ind[i];
        res.push_back( X[j] );
    }
    return res;
}


std::vector<dType> condsubset(std::vector<dType>& X, const std::vector<bool>& cond) {
    std::vector<dType> res;
    for (int i = 0; i < cond.size(); i++){
        if (cond[i])
            res.push_back(X[i]);
    }
    return res;
}

std::vector<iType> condsubset(std::vector<iType>& X, const std::vector<bool>& cond) {
    std::vector<iType> res;
    for (int i = 0; i < cond.size(); i++){
        if (cond[i])
            res.push_back(X[i]);
    }
    return res;
}


std::vector<dType> sum(const std::vector<std::vector<dType>>& x) {
    std::vector<dType> res(x.size());
    for (int i = 0; i < x.size(); i++) {
        res[i] = std::accumulate(x[i].begin(), x[i].end(), 0.0);
    }
    return res;
}


dType sum(const std::vector<dType>& x) {
    return std::accumulate(x.begin(), x.end(), 0.0);
}

dType sum(const std::vector<unsigned>& x) {
    return std::accumulate(x.begin(), x.end(), 0);
}

dType sum(const std::vector<dType>& x, const dType& val) {
    return std::accumulate(x.begin(), x.end(), val);;
}

dType sum(const std::vector<dType>& X, const std::vector<int>& ind) {
    dType res = 0;
    int j, i;
    #pragma omp parallel for private(i, j) reduction(+:res)
    for (i = 0; i < ind.size(); i++){
        j = ind[i];
        res += X[j];
    }
    return res;
}


std::vector<dType> pow(const std::vector<dType>& X, const dType& p) {
    std::vector<dType> Out;
    Out.reserve(X.size());
	for (int i = 0; i < X.size(); i++){
        Out.push_back( pow(X[i], p) );
    }
    return Out;
}

std::vector<dType> pow(const dType& X, const std::vector<dType>& p) {
    std::vector<dType> Out;
    Out.reserve(p.size());
    for (int i = 0; i < p.size(); i++){
        Out.push_back( pow(X, p[i]) );
    }
    return Out;
}

void Clear(std::vector<std::vector<dType>>& X) {
	for (int i = 0; i < X.size(); i++){
        X[i].clear();
    }
}

    //std::transform(myv1.begin(), myv1.end(), myv1.begin(), std::bind1st(std::multiplies<T>(),3));
    //std::transform(Array.begin(), Array.end(), Array.begin(), std::bind2nd(std::multiplies<dType>(), 0.5));
    //transform(v.begin(), v.end(), v.begin(), _1 * 3);
    //dType sum = std::accumulate(Array.begin(), Array.end(), 0.0);
    //std::transform(Array.begin(), Array.end(), Array.begin(), static_cast<dType (*)(dType)>(std::sqrt));
    //std::transform(a.begin(), a.end(), b.begin(), r.begin(), std::multiplies<float>());
