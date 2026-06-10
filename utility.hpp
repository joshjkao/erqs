#include <chrono>
#include <functional>
#include <pthread.h>
#include <type_traits>

template <typename ResultType> struct TimedResult {
  long time;
  ResultType ret;
};
template <> struct TimedResult<void> {
  long time;
};
auto time_in_ms(auto &&func) {
  using ResultType = std::invoke_result_t<decltype(func)>;
  auto start = std::chrono::high_resolution_clock::now();
  if constexpr (std::is_void_v<ResultType>) {
    std::invoke(std::forward<decltype(func)>(func));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return TimedResult<ResultType>{.time = duration.count()};
  } else {
    auto ret = std::invoke(std::forward<decltype(func)>(func));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return TimedResult<ResultType>{.time = duration.count(),
                                   .ret = std::move(ret)};
  }
}
