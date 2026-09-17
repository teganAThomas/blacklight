// Blacklight radiation integrator - unpolarized radiation integration

// C++ headers
#include <algorithm>  // min
#include <cmath>      // exp, expm1
#include <limits>     // numeric_limits

// Library headers
#include <omp.h>  // pragmas

// Blacklight headers
#include "radiation_integrator.hpp"
#include "../blacklight.hpp"         // Physics, enums
#include "../utils/array.hpp"        // Array

//--------------------------------------------------------------------------------------------------

// Function for integrating unpolarized radiative transfer equation
// Inputs: (none)
// Outputs: (none)
// Notes:
//   Assumes sample_num[adaptive_level], sample_len[adaptive_level], j_i[adaptive_level],
//       alpha_i[adaptive_level], and momentum_factors[adaptive_level] have been set.
//   Assumes sample_pos[adaptive_level] has been set if image_time == true or image_length == true.
//   Assumes sample_dir[adaptive_level] has been set if image_length == true.
//   Assumes cell_values[adaptive_level] has been set if image_lambda_ave == true or
//       image_emission_ave == true or image_tau_int == true.
//   Allocates and initializes image[adaptive_level].
//   Deallocates j_i[adaptive_level] and alpha_i[adaptive_level] if adaptive_level > 0.
//   Deallocates cell_values[adaptive_level] if render_num_images <= 0 and adaptive_level > 0.
void RadiationIntegrator::IntegrateUnpolarizedRadiation()
{
  // Allocate image array
  int num_pix = camera_num_pix;
  if (adaptive_level > 0)
    num_pix = block_counts[adaptive_level] * block_num_pix;
  if (first_time or adaptive_level > 0)
    image[adaptive_level].Allocate(image_num_quantities, num_pix);
  image[adaptive_level].Zero();

  // Calculate units
  double x_unit = Physics::gg_msun * mass_msun / (Physics::c * Physics::c);
  double t_unit = x_unit / Physics::c;

  // Work in parallel
  #pragma omp parallel
  {
    // Allocate scratch space
    double gcov[4][4];
    double gcon[4][4];

    // Go through frequencies and pixels
    #pragma omp for schedule(static) collapse(2)
    for (int l = 0; l < image_num_frequencies; l++)
      for (int m = 0; m < num_pix; m++)
      {
        // Extract number of steps
        int num_steps = sample_num[adaptive_level](m);

        // Prepare integrated quantities
        double integrated_lambda = 0.0;
        double integrated_emission = 0.0;
        double integrated_mcscat = 0.0;
        double x1_init = sample_pos[adaptive_level](m,0,1);
        double x2_init = sample_pos[adaptive_level](m,0,2);
        double x3_init = sample_pos[adaptive_level](m,0,3);
        bool plane_sign =
            camera_x[1] * x1_init + camera_x[2] * x2_init + camera_x[3] * x3_init > 0.0;
        int crossings_count = 0;

        double previous_delta_tau = 0.;
        int photosphere_steps = 0;
        // Go through samples
        for (int n = 0; n < num_steps; n++)
        {
          // Extract and calculate useful values
          double delta_lambda = sample_len[adaptive_level](m,n);
          double delta_lambda_cgs =
              delta_lambda * x_unit / (image_frequencies(l) * momentum_factors[adaptive_level](m));
          double t_cgs = sample_pos[adaptive_level](m,n,0) * t_unit;
          double x1 = sample_pos[adaptive_level](m,n,1);
          double x2 = sample_pos[adaptive_level](m,n,2);
          double x3 = sample_pos[adaptive_level](m,n,3);
          double kcov[4];
          kcov[0] = sample_dir[adaptive_level](m,n,0);
          kcov[1] = sample_dir[adaptive_level](m,n,1);
          kcov[2] = sample_dir[adaptive_level](m,n,2);
          kcov[3] = sample_dir[adaptive_level](m,n,3);
          double j = std::numeric_limits<double>::quiet_NaN();
          if (image_light or image_emission or image_emission_ave)
            j = j_i[adaptive_level](l,m,n);
          double alpha = std::numeric_limits<double>::quiet_NaN();
          if (image_light or image_tau or image_tau_int)
            alpha = alpha_i[adaptive_level](l,m,n);
          double ss = j / alpha;
          double delta_tau = alpha * delta_lambda_cgs;
          double exp_neg = std::exp(-delta_tau);
          double expm1 = std::expm1(delta_tau);
          bool optically_thin = delta_tau <= delta_tau_max;
          double mcscat;
          if(mc_error)
            mcscat = sample_scattering[adaptive_level](m,n,l);

          int trace_pixel = m;

          // Integrate light
          if (image_light)
          {
            if (alpha > 0.0)
            {
              if (optically_thin)
              {
                image[adaptive_level](l,m) = exp_neg * (image[adaptive_level](l,m) + ss * expm1);
                if(mc_error)
                  image[adaptive_level](image_offset_scat_err+l,m) = std::sqrt(std::pow(image[adaptive_level](image_offset_scat_err+l,m)*exp_neg,2.)+std::pow(exp_neg*expm1/alpha,2.)*scat_err[adaptive_level](l,m,n));
              }
              else 
              {
                image[adaptive_level](l,m) = ss;
                if(mc_error)
                  image[adaptive_level](image_offset_scat_err+l,m) = std::sqrt(scat_err[adaptive_level](l,m,n)/(alpha*alpha));
              }
            }
            else
            {
              image[adaptive_level](l,m) += j * delta_lambda_cgs;
              if(mc_error)
                image[adaptive_level](image_offset_scat_err+l,m) = std::sqrt(
                    std::pow(image[adaptive_level](image_offset_scat_err+l,m), 2.)
                    + std::pow(delta_lambda_cgs, 2.) * scat_err[adaptive_level](l,m,n));
            }
          }


          // Integrate alternative image quantities
          if (image_time and l == 0)
            image[adaptive_level](image_offset_time,m) =
                std::min(image[adaptive_level](image_offset_time,m), t_cgs);
          if (image_length and l == 0)
          {
            CovariantGeodesicMetric(x1, x2, x3, gcov);
            ContravariantGeodesicMetric(x1, x2, x3, gcon);
            double temp_a[4] = {};
            for (int a = 1; a < 4; a++)
              for (int mu = 0; mu < 4; mu++)
                temp_a[a] += (gcon[a][mu] - gcon[0][a] * gcon[0][mu] / gcon[0][0]) * kcov[mu];
            double dl_dlambda_sq = 0.0;
            for (int a = 1; a < 4; a++)
              for (int b = 1; b < 4; b++)
                dl_dlambda_sq += gcov[a][b] * temp_a[a] * temp_a[b];
            image[adaptive_level](image_offset_length,m) +=
                std::sqrt(dl_dlambda_sq) * delta_lambda * x_unit;
          }
          if (image_lambda or image_lambda_ave)
            integrated_lambda += delta_lambda_cgs;
          if (image_emission or image_emission_ave)
            integrated_emission += j * delta_lambda_cgs;
          if (image_tau)
            image[adaptive_level](image_offset_tau+l,m) += delta_tau;
          if (image_mcscat)
            integrated_mcscat += mcscat * delta_lambda_cgs;
          //the problem with how i do things is that mcscat is in comoving frame and 
          // tau is calculated in invariant form using alpha_I. 
          // dont know if we can simply convert alpha_I to alpha_nu and be fine and accurately integrate without changing affine parameter
          // I think it's fine because if it's invariant then it applies to any frame
          // the big thing is getting the proper comoving nu_cgs. for now i'm just using the camera frame frequency
          // i think the delta_lambda_cgs and its image_frequencies in the denom already handles converting alpha_I potentially
          if (image_mcscat_ave){
            double tau_0 = 0.3;
            double tau_1 = 3.0;

            if(image[adaptive_level](image_offset_tau+l,m) >= tau_0 and (image[adaptive_level](image_offset_tau+l,m) <= tau_1 || previous_delta_tau < tau_0) and delta_tau!=0.0){
              image[adaptive_level](image_offset_mcscat_ave+l,m) += mcscat * delta_tau;
              image[adaptive_level](image_offset_mcscat_ave_weight+l,m) += delta_tau;
            }
          }


          if (image_lambda_ave and not std::isnan(cell_values[adaptive_level](0,m,n)))
            for (int a = 0; a < CellValues::num_cell_values; a++)
            {
              int index = image_offset_lambda_ave + l * CellValues::num_cell_values + a;
              image[adaptive_level](index,m) +=
                  cell_values[adaptive_level](a,m,n) * delta_lambda_cgs;
            }
          if (image_emission_ave and not std::isnan(cell_values[adaptive_level](0,m,n)))
            for (int a = 0; a < CellValues::num_cell_values; a++)
            {
              int index = image_offset_emission_ave + l * CellValues::num_cell_values + a;
              image[adaptive_level](index,m) +=
                  cell_values[adaptive_level](a,m,n) * j * delta_lambda_cgs;
            }
          if (image_tau_int and not std::isnan(cell_values[adaptive_level](0,m,n)))
          {
            if (optically_thin)
              for (int a = 0; a < CellValues::num_cell_values; a++)
              {
                int index = image_offset_tau_int + l * CellValues::num_cell_values + a;
                image[adaptive_level](index,m) = exp_neg
                    * (image[adaptive_level](index,m) + cell_values[adaptive_level](a,m,n) * expm1);
              }
            else
              for (int a = 0; a < CellValues::num_cell_values; a++)
              {
                int index = image_offset_tau_int + l * CellValues::num_cell_values + a;
                image[adaptive_level](index,m) = cell_values[adaptive_level](a,m,n);
              }
          }
          if (image_photosphere_int and not std::isnan(cell_values[adaptive_level](0,m,n)))
          {
            double tau_0 = 0.5;
            double tau_1 = 1.7;
            for (int a = 0; a < CellValues::num_cell_values; a++){
              //int index = image_offset_photosphere_int + l * (CellValues::num_cell_values+3) + a;
              int index = image_offset_photosphere_int + l * (CellValues::num_cell_values) + a;
              //delta_tau sometimes equals 0 which causes Nan's
              if(image[adaptive_level](image_offset_tau+l,m) >= tau_0 and (image[adaptive_level](image_offset_tau+l,m) <= tau_1 || previous_delta_tau < tau_0) and (not std::isnan(cell_values[adaptive_level](a,m,n))) and delta_tau!=0.0){
                image[adaptive_level](index,m)+= cell_values[adaptive_level](a,m,n);///delta_tau;
              }
            }
            /*int index = image_offset_photosphere_int + l * (CellValues::num_cell_values+3) +CellValues::num_cell_values;
            if(image[adaptive_level](image_offset_tau+l,m) >= tau_0 and (image[adaptive_level](image_offset_tau+l,m) <= tau_1 || previous_delta_tau < tau_0) and (not std::isnan(x1)) and delta_tau!=0.0){
                image[adaptive_level](index,m)+= x1;///delta_tau;
              }else{
                image[adaptive_level](index,m) += 0.0;
              }
            index++;
            if(image[adaptive_level](image_offset_tau+l,m) >= tau_0 and (image[adaptive_level](image_offset_tau+l,m) <= tau_1 || previous_delta_tau < tau_0) and (not std::isnan(x2)) and delta_tau!=0.0){
                image[adaptive_level](index,m)+= x2;///delta_tau;
              }else{
                image[adaptive_level](index,m) += 0.0;
              }
            index++;
            if(image[adaptive_level](image_offset_tau+l,m) >= tau_0 and (image[adaptive_level](image_offset_tau+l,m) <= tau_1 || previous_delta_tau < tau_0) and (not std::isnan(x3)) and delta_tau!=0.0){
                image[adaptive_level](index,m)+= x3;///delta_tau;
              }else{
                image[adaptive_level](index,m) += 0.0;
              }*/
            //image[adaptive_level](index,m)+= x1;
            //image[adaptive_level](index+1,m)+= x2;
            //image[adaptive_level](index+2,m)+= x3;
            photosphere_steps++;
            
          }
          if (image_crossings and l == 0)
          {
            bool plane_sign_new = camera_x[1] * x1 + camera_x[2] * x2 + camera_x[3] * x3 > 0.0;
            if (plane_sign_new != plane_sign)
              crossings_count++;
            plane_sign = plane_sign_new;
          }
          if(image_photosphere_int || image_mcscat_ave){
            previous_delta_tau = image[adaptive_level](image_offset_tau+l,m);
          }
        }
        
        // Store integrated quantities
        if (image_lambda)
          image[adaptive_level](image_offset_lambda+l,m) = integrated_lambda;
        if (image_emission)
          image[adaptive_level](image_offset_emission+l,m) = integrated_emission;
        if (image_mcscat)
          image[adaptive_level](image_offset_mcscat+l,m) = integrated_mcscat;
        if (image_crossings and l == 0)
          image[adaptive_level](image_offset_crossings,m) = static_cast<double>(crossings_count);

        // Normalize integrated quantities
        if (image_mcscat_ave){
          image[adaptive_level](image_offset_mcscat_ave+l,m) /= image[adaptive_level](image_offset_mcscat_ave_weight+l,m);
        }
        if (image_lambda_ave)
          for (int a = 0; a < CellValues::num_cell_values; a++)
          {
            int index = image_offset_lambda_ave + l * CellValues::num_cell_values + a;
            image[adaptive_level](index,m) /= integrated_lambda;
          }
        if (image_emission_ave)
          for (int a = 0; a < CellValues::num_cell_values; a++)
          {
            int index = image_offset_emission_ave + l * CellValues::num_cell_values + a;
            image[adaptive_level](index,m) /= integrated_emission;
          }
        if (image_photosphere_int){
          for (int a = 0; a < CellValues::num_cell_values; a++)
          {
            //int index = image_offset_photosphere_int + l * (CellValues::num_cell_values+3) + a;
            int index = image_offset_photosphere_int + l * (CellValues::num_cell_values) + a;
            image[adaptive_level](index,m) /= photosphere_steps;
          }
          /*int index = image_offset_photosphere_int + l * (CellValues::num_cell_values+3)+CellValues::num_cell_values;
          image[adaptive_level](index,m) /= photosphere_steps;
          index++;
          image[adaptive_level](index,m) /= photosphere_steps;
          index++;
          image[adaptive_level](index,m) /= photosphere_steps;*/
        }   
    }
    // Transform I_nu/nu^3 to I_nu
    if (image_light)
    {
      #pragma omp for schedule(static) collapse(2)
      for (int l = 0; l < image_num_frequencies; l++)
        for (int m = 0; m < num_pix; m++)
        {
          double nu_cu = image_frequencies(l) * image_frequencies(l) * image_frequencies(l);
          image[adaptive_level](l,m) *= nu_cu;
          if(mc_error)
            image[adaptive_level](image_offset_scat_err+l,m) *= nu_cu;
        }
    }
  }

  // Free memory
  if (adaptive_level > 0)
  {
    j_i[adaptive_level].Deallocate();
    alpha_i[adaptive_level].Deallocate();
    if (render_num_images <= 0)
      cell_values[adaptive_level].Deallocate();
    if(mc_error)
      scat_err[adaptive_level].Deallocate();
  }
  return;
}
