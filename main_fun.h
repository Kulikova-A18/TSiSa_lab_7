#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <vector>
#include <random>
#include <algorithm>

#include "main_graph.h"
/*
    Структура для описания условий задачи
*/
struct Data
{
    std::pair<double,double> X_interval = {0.0, 3.141592653589793238463}; // Interval
    size_t K_ = 100; // Number of samples
    double amplitude_of_noise = 0.5; // Uniform noise amplitude
    double P = 0.95; // Probability of getting into the vicinity of an extremum
    double eps = 0.01; // Uncertainty interval
    double J = 90000000.0;
    size_t L_ = 10u; // Uncertainty interval
    std::vector <size_t> radius_ {3,5}; // Sliding window size
    std::vector <size_t> M_;
    std::vector <double> alpha;
    std::vector <double> xk;
    std::vector <double> F_xk_;
    std::vector <double> F_xk_noised;
    std::vector <double> F_xk_filtered;
    std::vector <double> lambda; // Convolution weights
    std::pair<double,double> omega_delta{90000000.0,90000000.0};

    void init() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution <double> dis (-amplitude_of_noise/2.0,amplitude_of_noise/2.0);

        for (size_t index = 0u;index <= L_;++index)
            lambda.push_back(0.1 * index);
        // Calculate xk
        for (size_t index = 0u;index <= K_;++index)
            xk.push_back(calculate_x_k(index));
        // Calculate f(xk) and f(xk) with noise
        for (size_t index = 0u;index <= K_;++index) {
            F_xk_.push_back(signal(xk[index]));
            F_xk_noised.push_back(signal_plus_noise(xk[index],dis(gen)));
        }
        // Calculate M
        for (size_t index = 0u;index < radius_.size(); ++index)
            M_.push_back((radius_[index]-1u)/2u);

        std::cout<<"";
    }

    const double calculate_x_k(size_t k) {
        if (k >= 0u && k <= 100u)
            return X_interval.first + k * (X_interval.second - X_interval.first)/K_;
        else
            throw std::runtime_error("Неверное значение для k!");
    }
    /**
        Функция вычисления сигнала для xk
    */
    const double signal(double Xk_) {
        return std::sin(Xk_) + 0.5;
    }
    /**
        Функция вычисления сигнала с шумом для xk
    */
    const double signal_plus_noise(double Xk_,double noise) {
        return std::sin(Xk_) + 0.5 + noise;
    }
    /**
        Функция получения числа испытаний N
    */
    const double get_N()
        {
            return (log(1.0-P)/log(1.0 - ( eps/( X_interval.second - X_interval.first ) )));
        }
    /**
        Функция вычисления среднего геометрического
    */
    double geom_average(std::vector<double> alpha,size_t M,size_t k) {
        double temp = 1.0;
        if(k < M || k > K_ - M)	// k-M должно быть больше или равно 0
            {
                return 0.0;
            }
        size_t degree = 0u;

        for (size_t index = k - M ; index <= k + M ;++index)
            {
                degree = index + M - k ;
                if (degree >= alpha.size() )
                    {
                        degree  = alpha[ alpha.size()-(degree-(alpha.size()-1u))-1u ];
                    }
                temp *= pow(F_xk_noised[index],alpha[degree]);
            }
        return temp;
    }
    double sum(std::vector<double> alpha,size_t a,size_t b) {
        double temp = 0.0;
        for (size_t index = a;index <= b;++index)
            {
                if (index + 1u > alpha.size())
                    {
                        temp+=alpha[index-alpha.size()];
                        continue;
                    }
                temp+=alpha[index];
            }
        return temp;
    }
    void init_Fx_filtered(std::vector <double> input_alpha,size_t M) {
        F_xk_filtered.clear();
        for (size_t index = 0u;index <= K_;++index)
            {
                F_xk_filtered.push_back(geom_average(input_alpha,M,index));
            }
    }
};

class main_fun_worker {
    double N_;
    V10 v10;
    std::vector <Save> save;
    Data data;
    public:
        main_fun_worker() {
            data.init();
            N_ = data.get_N();
        }
        void find_best();
        void pass();

};
