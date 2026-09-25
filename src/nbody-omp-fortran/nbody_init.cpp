// SPDX-License-Identifier: CC0-1.0
#include <iostream>
#include <random>

extern "C" {

struct particle {
  float pos[3];
  float vel[3];
  float acc[3];
  float mass;
};

void nbody_init_particles(particle *particles, int n) {
  std::mt19937 gen_pos(42);
  std::uniform_real_distribution<float> unif_pos(0.0f, 1.0f);

  for (int i = 0; i < n; ++i) {
    particles[i].pos[0] = unif_pos(gen_pos);
    particles[i].pos[1] = unif_pos(gen_pos);
    particles[i].pos[2] = unif_pos(gen_pos);
  }

  std::mt19937 gen_vel(42);
  std::uniform_real_distribution<float> unif_vel(-1.0f, 1.0f);

  for (int i = 0; i < n; ++i) {
    particles[i].vel[0] = unif_vel(gen_vel) * 1.0e-3f;
    particles[i].vel[1] = unif_vel(gen_vel) * 1.0e-3f;
    particles[i].vel[2] = unif_vel(gen_vel) * 1.0e-3f;
  }

  for (int i = 0; i < n; ++i) {
    particles[i].acc[0] = 0.0f;
    particles[i].acc[1] = 0.0f;
    particles[i].acc[2] = 0.0f;
  }

  std::mt19937 gen_mass(42);
  std::uniform_real_distribution<float> unif_mass(0.0f, 1.0f);
  const float scale = static_cast<float>(n);

  for (int i = 0; i < n; ++i) {
    particles[i].mass = scale * unif_mass(gen_mass);
  }
}

void nbody_print_summary(float kenergy, double total_time, double av,
                         double dev) {
  std::cout << "\n";
  std::cout << "# Total Energy        : " << kenergy << "\n";
  std::cout << "# Total Time (s)      : " << total_time << "\n";
  std::cout << "# Average Performance : " << av << " +- " << dev << "\n";
  std::cout << "===============================\n" << std::flush;
}

}
