#include "main_fun.h"

void main_fun_worker::pass()
{
    Save temp_save;
    std::pair<double,double> temp_omega_delta;
    double temp_J = 0.0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution <double> dis (0.0, 1.0);
    std::vector<double> temp_alpha;

    double temp;

    temp_save.graphic.x = data.xk;
    temp_save.graphic.Function = data.F_xk_;
    temp_save.graphic.Function_noised = data.F_xk_noised;

    for (size_t jump = 0u; jump < data.radius_.size(); ++jump)
    {
        double r = data.radius_[jump];
        double M = data.M_[jump];

        temp_save.r = data.radius_[jump];
        // for different weight values
        for (size_t l = 0u; l < data.lambda.size(); ++l)
        {
            temp_save.h = data.lambda[l];
            // Number of tests N
            for (size_t index = 0u; index < N_; ++index)
            {
                temp_alpha.reserve(M+1u);
                temp_alpha.resize(M+1u);
                temp_alpha.back() = dis(gen);
                if (data.M_[jump] >= 2u)
                {
                    for (size_t count = 1u; count < data.M_[jump] ; ++count)
                    {
                        std::uniform_real_distribution <double> temp_dis(0.0,1.0 - data.sum(temp_alpha,M, r-M-1));
                        temp_alpha.at(count) = 0.5 * temp_dis(gen);
                    }
                }
                temp_alpha.front() = 0.5*(1.0 - data.sum(temp_alpha,1u,r-2u));
                //***Already generated alpha set ***

                // Find the filtered function
                data.init_Fx_filtered(temp_alpha,M);

                // ******
                // We find the noise criterion (according to “Manhattan”)
                double summ = 0.0;
                for (size_t kol = 1u; kol <= data.K_; ++kol)
                {
                    summ += std::abs(data.F_xk_filtered[kol] - data.F_xk_filtered[kol - 1u]);
                }
                temp_omega_delta.first = sqrt(summ); summ=0.0;
                // Find the difference criterion (according to “Manhattan”)
                for (size_t kol = 0u; kol <= data.K_; ++kol)
                {
                    summ += std::abs(data.F_xk_filtered[kol] - data.F_xk_noised[kol]);
                }
                temp_omega_delta.second = sqrt(summ/data.K_);

                temp_J = ( data.lambda[l] * temp_omega_delta.first) + (1 - data.lambda[l])*temp_omega_delta.second;

                if (temp_J < data.J)
                {
                    data.J = temp_J;
                    data.omega_delta =	temp_omega_delta;
                    temp_save.alpha = temp_alpha;
                    temp_save.metrics.delta = temp_omega_delta.second;
                    temp_save.metrics.omega = temp_omega_delta.first;
                    temp_save.metrics.J = temp_J;
                    temp_save.graphic.Function_filtered = data.F_xk_filtered;
                    // |delta| + |omega|
                    temp_save.distance = std::abs(temp_save.metrics.delta) + std::abs(temp_save.metrics.omega);
                }
            }
            data.J = 90000000.0;
            data.omega_delta = {90000000.0,90000000.0};
            save.push_back(temp_save);
        }

        if (jump == 0u)
        {
            v10.resut_r1.saves = save ;
            save.clear();
        }
        else
        {
            v10.resut_r2.saves = save ;
            save.clear();
        }
    }

    find_best();
    std::cout<<"\t ************************************ \n";
    std::cout<<"\t ***** FOR SLIDING WINDOW R = 3 ***** \n";
    std::cout<<"\t ************************************ \n\n";
    v10.resut_r1.print();

    std::cout<<"\t ************************************ \n";
    std::cout<<"\t ***** FOR SLIDING WINDOW R = 5 ***** \n";
    std::cout<<"\t ************************************ \n\n";
    v10.resut_r2.print();
}
void main_fun_worker::find_best()
{
    Save temp_best;
    double best_distance = 999999900.0;
    for (Save save:v10.resut_r1.saves)
    {
        if (save.distance < best_distance)
            {
                best_distance = save.distance;
                temp_best = save;
            }
    }
    v10.resut_r1.best_save = temp_best;
    best_distance = 999999900.0;
    for (Save save:v10.resut_r2.saves)
    {
        if (save.distance < best_distance)
            {
                best_distance = save.distance;
                temp_best = save;
            }
    }
    v10.resut_r2.best_save = temp_best;
}
