#pragma once

#include <vector>
#include <unordered_map>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <future>
#include <thread>
#include <atomic>
#include <utility>
#include "gsi/CLL.hpp"

namespace fmx::gsi {

struct CLLQuantization {
  double theta_deg_step{2.0};
  double tau_step{0.1};
  double alpha_step{0.25};
  std::vector<double> Ma_anchors{0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0};
  int sample_count{4000};
  int workers{1};
  int gh_order{8};
};

class CLLRuntime {
public:
  explicit CLLRuntime(const CLLQuantization& q = {});
  ~CLLRuntime();

  // Non-copyable
  CLLRuntime(const CLLRuntime&) = delete;
  CLLRuntime& operator=(const CLLRuntime&) = delete;

  std::pair<double,double> query(double theta, double Ma, double tau, double an, double at) const;

private:
  struct Key {
    double theta, Ma, tau, an, at;
    bool operator==(const Key& o) const {
      return theta==o.theta && Ma==o.Ma && tau==o.tau && an==o.an && at==o.at;
    }
  };

  struct KeyHash {
    size_t operator()(const Key& k) const {
      auto h1 = std::hash<long long>{}(static_cast<long long>(k.theta*1e6));
      auto h2 = std::hash<long long>{}(static_cast<long long>(k.Ma*1e6));
      auto h3 = std::hash<long long>{}(static_cast<long long>(k.tau*1e6));
      auto h4 = std::hash<long long>{}(static_cast<long long>(k.an*1e6));
      auto h5 = std::hash<long long>{}(static_cast<long long>(k.at*1e6));
      return (((h1*1315423911u) ^ h2)*2654435761u ^ h3) ^ (h4*1099511628211ull) ^ h5;
    }
  };

  using Promise = std::promise<std::pair<double,double>>;
  using Future = std::shared_future<std::pair<double,double>>;

  struct Entry { Future fut; };

  CLLQuantization q_;
  mutable std::unordered_map<Key, Entry, KeyHash> cache_;
  mutable std::mutex mtx_;

  mutable std::deque<std::pair<Key, std::shared_ptr<Promise>>> queue_;
  mutable std::mutex q_mtx_;
  mutable std::condition_variable q_cv_;
  std::atomic<bool> stop_{false};
  std::vector<std::thread> workers_;

  Key quantize(double theta, double Ma, double tau, double an, double at) const;
  void worker_loop();
  static std::pair<double,double> compute_cll(double theta, double Ma, double tau, double an, double at, int gh_order);
};

} // namespace fmx::gsi
