#ifndef RANDOM_THIJS
#define RANDOM_THIJS
#include <random>
#include <chrono>
#include <thread>
#include <type_traits>

struct rnd_t {
  std::mt19937 rndgen;

  rnd_t() {
    std::mt19937 rndgen_t(get_seed());
    rndgen = rndgen_t;
  }

  rnd_t(size_t seed) {
    std::mt19937 rndgen_t(seed);
    rndgen = rndgen_t;
  }

  int get_seed() {
    const auto tt = static_cast<int64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    //auto tid = std::this_thread::get_id();
    //const uint64_t e3{ std::hash<std::remove_const_t<decltype(tid)>>()(tid) };
    const auto e3 = 0.0;
    auto output = static_cast<int>(tt + e3);
    if (output < 0) output *= -1;
    return output;
  }

  std::uniform_real_distribution<float> unif_dist =
    std::uniform_real_distribution<float>(0.0f, 1.0f);


  std::uniform_int_distribution<> world_dist;

  size_t random_number(size_t n)    {
    if(n <= 1) return 0;
    return static_cast<size_t>(std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(rndgen));
  }

  size_t random_pos() {
    return static_cast<size_t>(world_dist(rndgen));
  }

  void set_world_size(size_t world_size) {
    world_dist = std::uniform_int_distribution<>(0, static_cast<int>(world_size) - 1);
  }

  float uniform()    {
    return unif_dist(rndgen);
  }

  void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
  }

  bool bernouilli(double p) {
    std::bernoulli_distribution d(p);
    return(d(rndgen));
  }
};


#endif /* rand_t.h */
