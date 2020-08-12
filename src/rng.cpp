#include "cglass/rng.hpp"

RNG::RNG(unsigned long seed) {
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  rng_ = gsl_rng_alloc(T);
  gsl_rng_set(rng_, seed);
}

RNG::~RNG() { gsl_rng_free(rng_); }

unsigned long RNG::GetSeed() const { return gsl_rng_get(rng_); }

const long RNG::RandomInt(const long n) { return gsl_rng_uniform_int(rng_, n); }

void *RNG::GetState() { return gsl_rng_state(rng_); }

size_t RNG::GetSize() { return gsl_rng_size(rng_); }

const double RNG::RandomUniform() { return gsl_rng_uniform_pos(rng_); }

const int RNG::RandomPoisson(const double mean) {
  return gsl_ran_poisson(rng_, mean);
}

const double RNG::RandomNormal(const double sigma) {
  return gsl_ran_gaussian_ziggurat(rng_, sigma);
}

void RNG::RandomUnitVector(const int n_dim, double *vec) {
  double w = 1.0;
  if (n_dim == 3) {
    double z = 2.0 * gsl_rng_uniform_pos(rng_) - 1.0;
    w = sqrt(1 - z * z);
    vec[2] = z;
  }

  double t = 2.0 * M_PI * gsl_rng_uniform_pos(rng_);
  double x = w * cos(t);
  double y = w * sin(t);
  vec[0] = x;
  vec[1] = y;
  Logger::Trace("Generated random unit vector: [%2.2f %2.2f %2.2f]", vec[0],
                vec[1], vec[2]);
}
void RNG::RandomCoordinate(const space_struct *const s, double *vec,
                           const double buffer) {
  double R = s->radius;
  int n_dim = s->n_dim;
  if (R - buffer < 0) {
    Logger::Error(
        "RNG tried to generate a random coordinate but the buffer "
        "received %2.2f is larger than the system radius %2.2f",
        R, buffer);
  }
  double mag;
  switch (s->type) {
    // If no boundary, insert wherever
    case +boundary_type::none:  // none
      for (int i = 0; i < n_dim; ++i) {
        vec[i] = (2.0 * gsl_rng_uniform_pos(rng_) - 1.0) * (R - buffer);
      }
      break;
    // box type boundary
    case +boundary_type::box:  // box
      for (int i = 0; i < n_dim; ++i) {
        vec[i] = (2.0 * gsl_rng_uniform_pos(rng_) - 1.0) * (R - buffer);
      }
      break;
    // spherical boundary
    case +boundary_type::sphere:  // sphere
      RandomUnitVector(n_dim, vec);
      mag = gsl_rng_uniform_pos(rng_) * (R - buffer);
      for (int i = 0; i < n_dim; ++i) {
        vec[i] *= mag;
      }
      break;
    // budding yeast boundary type
    case +boundary_type::budding:  // budding
    {
      double r = s->bud_radius;
      double roll = gsl_rng_uniform_pos(rng_);
      double v_ratio = 0;
      if (n_dim == 2) {
        v_ratio = SQR(r) / (SQR(r) + SQR(R));
      } else {
        v_ratio = CUBE(r) / (CUBE(r) + CUBE(R));
      }
      mag = gsl_rng_uniform_pos(rng_);
      RandomUnitVector(n_dim, vec);
      if (roll < v_ratio) {
        // Place coordinate in daughter cell
        mag *= (r - buffer);
        for (int i = 0; i < n_dim; ++i) {
          vec[i] *= mag;
        }
        vec[n_dim - 1] += s->bud_height;
      } else {
        mag *= (R - buffer);
        for (int i = 0; i < n_dim; ++i) {
          vec[i] *= mag;
        }
      }
      break;
    }
    // For wall, insert wherever
    case +boundary_type::wall:
    {
      for (int i = 0; i < n_dim; ++i) {
        vec[i] = (2.0 * gsl_rng_uniform_pos(rng_) - 1.0) * (R - buffer);
      }
      break;
    }
    default:
      Logger::Error("Boundary type unrecognized in RandomCoordinate");
  }
  Logger::Trace("Generated random coordinate: [%2.2f %2.2f %2.2f]", vec[0],
                vec[1], vec[2]);
}

void RNG::RandomBoundaryCoordinate(const space_struct *const s, double *vec) {
  double R = s->radius;
  int n_dim = s->n_dim;
  switch (s->type) {
    case +boundary_type::none: {
      Logger::Error(
          "Random boundary vector cannot be created for boundary type = none");
      break;
    }
    case +boundary_type::box: {
      int roll_pm = RandomInt(2) * 2 - 1;
      int roll_plane = RandomInt(n_dim);
      for (int i = 0; i < n_dim; ++i) {
        if (i == roll_plane) {
          vec[i] = roll_pm * R; // set one coordinate to +/- R
        } else {
          vec[i] = (2.0 * RandomUniform() - 1.0) * R; // others random
        }
      }
      break;
    }
    case +boundary_type::sphere: {
      RandomUnitVector(n_dim, vec);
      for (int i = 0; i < n_dim; ++i) {
        vec[i] *= R;
      }
      break;
    }
    case +boundary_type::budding: {
      double r = s->bud_radius;
      double d = s->bud_height;
      
      // alpha = angle of intersect on mother cell
      double alpha = acos((SQR(d) - SQR(r) + SQR(R)) / (2.0 * d * R));
      
      // beta = angle of intersect on daughter cell 
      double beta = acos((SQR(d) + SQR(r) - SQR(R)) / (2.0 * d * r)); 
      double a_ratio = 0;
      double roll = RandomUniform();
      if (n_dim == 2) {
        // arc length for segment  = 2*r*(pi-theta)
        a_ratio = r * (M_PI - beta) / (r * (M_PI - beta) + R * (M_PI - alpha));
      }
      else {
        // surface area for segment  = 2*pi*r^2*(1+cos(theta))
        a_ratio = 
            SQR(r) * (1 + cos(beta)) / (SQR(r) * (1 + cos(beta)) + SQR(R) * (1 + 
            cos(alpha)));
      }
      std::fill(vec, vec + n_dim, 0.0);
      double theta = 0.0;
      double rho = 1.0;
      if (roll < a_ratio) {
        // place on daughter cell boundary
        theta = RandomUniform() * (M_PI - beta); // 0 to pi-beta
        rho = r;
        vec[n_dim - 1] += d;
      } else {
        theta = RandomUniform() * (M_PI - alpha) + alpha; // alpha to pi
        rho = R;        
      }
      if (n_dim == 2) {
        vec[0] += rho * sin(theta) * (RandomInt(2) * 2 - 1); // choose +/- x side
        vec[1] += rho * cos(theta);
      } else {
        double phi = RandomUniform() * 2 * M_PI;
        vec[0] += rho * sin(theta) * cos(phi);
        vec[1] += rho * sin(theta) * sin(phi);
        vec[2] += rho * cos(theta);
      }
      break;
    }
    default: 
      Logger::Error("Boundary type unrecognized in RandomBoundaryCoordinate");
  }
  Logger::Trace("Generated random coordinate: [%2.2f %2.2f %2.2f]", vec[0],
                vec[1], vec[2]);
}
